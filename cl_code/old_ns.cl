#define radius_segments 30

int searchForNeighbors(
		int searchCell_,
		__global int * b_grid_cell_id_list,
		float4 position_,
		int myParticleId,
		__global struct particle * particles,
		__global struct extend_particle * ext_particles,
		int spaceLeft,
		float h,
		float simulationScale,
		int mode,
		int * radius_distrib,
		float r_thr
		)
{
	int baseParticleId = b_grid_cell_id_list[ searchCell_ ];
	int foundCount = 0;
	bool loop = true;
	int i = 0,j;
	float _distance,_distanceSquared;
	float r_thr_Squared = r_thr*r_thr;
	float2 neighbor_data;
	int myOffset;
	if(baseParticleId == -1){
		printf("\nBAD ASS\n");
		return foundCount;
	}
	int neighborParticleId = baseParticleId;
	if(spaceLeft>0)
		while( particles[neighborParticleId].cell_id == searchCell_){
			if (myParticleId != neighborParticleId)
			{
				float4 d = position_ - particles[neighborParticleId].pos;
				d.w = 0.0f;
				_distanceSquared = d.x * d.x + d.y * d.y + d.z * d.z; // inlined openCL dot(d,d)

				if( _distanceSquared <= r_thr_Squared )
				{
					_distance = SQRT( _distanceSquared );
					j = (int)(_distance*radius_segments/h);
					if(j<radius_segments) radius_distrib[j]++;

					// searchForNeighbors runs twice
					// first time with mode = 0, to build distribution
					// and 2nd time with mode = 1, to select 32 nearest neighbors
					if(mode)
					{
						myOffset = NEIGHBOUR_COUNT - spaceLeft + foundCount;
						if(myOffset >= NEIGHBOUR_COUNT) break;// New line fixing the bug with indeterminism. A. Palyanov 22.02.2013
						neighbor_data.x = neighborParticleId;
						neighbor_data.y = _distance * simulationScale; // scaled, OK
						ext_particles[myParticleId - 0/*OFFSET*/].neighbour_list[myOffset][0] = neighbor_data.x;
						ext_particles[myParticleId - 0/*OFFSET*/].neighbour_list[myOffset][1] = neighbor_data.y;
						//ext_particles[ myParticleId ]. + myOffset ] = neighbor_data;
						foundCount++;
					}
				}
			}
			++neighborParticleId;
		}
	return foundCount;
}




__kernel void k_neighbour_search(
		__global struct extend_particle * ext_particles,
		__global struct
				particle * particles,
        __global int * b_grid_cell_id_list,
		uint grid_cells_X,
		uint grid_cells_Y,
		uint grid_cells_Z,
		uint grid_cell_count,
        uint grid_offset,
		float h,
		float hashGridCellSize,
		float hashGridCellSizeInv,
		float simulationScale,
		float xmin,
		float ymin,
		float zmin,
		int PARTICLE_COUNT,
		int OFFSET,
		int LIMIT
){
	int id = get_global_id( 0 );
	if(id >= PARTICLE_COUNT || id + OFFSET > LIMIT){
		return;
	}
	id += OFFSET;
	float4 position_ = particles[ id ].pos;
	int myCellId = particles[id].cell_id;//& 0xffffff;// truncate to low 16 bits
	int searchCells[8];
	float r_thr2 = h * h;
	float closest_distances[NEIGHBOUR_COUNT];
	int closest_indexes[NEIGHBOUR_COUNT];
	int found_count = 0;
	for(int k=0;k<NEIGHBOUR_COUNT;++k){
		closest_distances[k] = r_thr2;
		closest_indexes[k] = -1;
	}
	searchCells[0] = myCellId - grid_offset;
	// p is the current particle position within the bounds of the hash grid
	float4 p;
	float4 p0 = (float4)( xmin, ymin, zmin, 0.0f );
	p = position_ - p0;

	// cf is the min,min,min corner of the current cell
	int4 cellFactors_ = cellFactors( position_, xmin, ymin, zmin, hashGridCellSizeInv );
	float4 cf;
	cf.x = cellFactors_.x * hashGridCellSize;
	cf.y = cellFactors_.y * hashGridCellSize;
	cf.z = cellFactors_.z * hashGridCellSize;

	// lo.A is true if the current position is in the low half of the cell for dimension A
	//float4 ttt = ( p - cf );
	int4 lo;
	int debug_id = 1455;
	lo = isless(p, h + cf);
	int p_id = id;//particles[id].particle_id;
	int4 delta;
    int4 one = (int4)( 1, 1, 1, 1 );
	delta = one + 2 * lo;
	int last_farthest = 0;

	searchCells[1] = searchCell( myCellId, delta.x, 0,       0,       grid_cells_X, grid_cells_Y, grid_cells_Z, grid_cell_count ) - grid_offset;
	searchCells[2] = searchCell( myCellId, 0,       delta.y, 0,       grid_cells_X, grid_cells_Y, grid_cells_Z, grid_cell_count ) - grid_offset;
	searchCells[3] = searchCell( myCellId, 0,       0,       delta.z, grid_cells_X, grid_cells_Y, grid_cells_Z, grid_cell_count ) - grid_offset;
	searchCells[4] = searchCell( myCellId, delta.x, delta.y, 0,       grid_cells_X, grid_cells_Y, grid_cells_Z, grid_cell_count ) - grid_offset;
	searchCells[5] = searchCell( myCellId, delta.x, 0,       delta.z, grid_cells_X, grid_cells_Y, grid_cells_Z, grid_cell_count ) - grid_offset;
	searchCells[6] = searchCell( myCellId, 0,       delta.y, delta.z, grid_cells_X, grid_cells_Y, grid_cells_Z, grid_cell_count ) - grid_offset;
	searchCells[7] = searchCell( myCellId, delta.x, delta.y, delta.z, grid_cells_X, grid_cells_Y, grid_cells_Z, grid_cell_count ) - grid_offset;

	int radius_distrib[radius_segments];
	int i=0,j;
	float r_thr = h;
	int mode = 0;
	int distrib_sum = 0;

	while( i<radius_segments )
	{
		radius_distrib[i]=0;
		i++;
	}

	while( mode<2 )
	{
		// search surrounding cell 1

		found_count += searchForNeighbors( searchCells[0], b_grid_cell_id_list, position_,
		                                  id, particles, ext_particles, NEIGHBOUR_COUNT - found_count,
		                                  h, simulationScale, mode, radius_distrib, r_thr  );

		found_count += searchForNeighbors( searchCells[1], b_grid_cell_id_list, position_,
		                                  id, particles, ext_particles, NEIGHBOUR_COUNT - found_count,
		                                  h, simulationScale, mode, radius_distrib, r_thr  );

		found_count += searchForNeighbors( searchCells[2], b_grid_cell_id_list, position_,
		                                  id, particles, ext_particles, NEIGHBOUR_COUNT - found_count,
		                                  h, simulationScale, mode, radius_distrib, r_thr  );

		found_count += searchForNeighbors( searchCells[3], b_grid_cell_id_list, position_,
		                                  id, particles, ext_particles, NEIGHBOUR_COUNT - found_count,
		                                  h, simulationScale, mode, radius_distrib, r_thr  );

		found_count += searchForNeighbors( searchCells[4], b_grid_cell_id_list, position_,
		                                  id, particles, ext_particles, NEIGHBOUR_COUNT - found_count,
		                                  h, simulationScale, mode, radius_distrib, r_thr  );

		found_count += searchForNeighbors( searchCells[5], b_grid_cell_id_list, position_,
		                                  id, particles, ext_particles, NEIGHBOUR_COUNT - found_count,
		                                  h, simulationScale, mode, radius_distrib, r_thr  );

		found_count += searchForNeighbors( searchCells[6], b_grid_cell_id_list, position_,
		                                  id, particles, ext_particles, NEIGHBOUR_COUNT - found_count,
		                                  h, simulationScale, mode, radius_distrib, r_thr );

		found_count += searchForNeighbors( searchCells[7], b_grid_cell_id_list, position_,
		                                  id, particles, ext_particles, NEIGHBOUR_COUNT - found_count,
		                                  h, simulationScale, mode, radius_distrib, r_thr );

		if(mode==0)
		{
			j=0;
			while(j<radius_segments)
			{
				distrib_sum += radius_distrib[j];
				if(distrib_sum==NEIGHBOUR_COUNT)
					break;
				if(distrib_sum> NEIGHBOUR_COUNT) {
					j--;
					break;
				}
				j++;
			}

			r_thr = (j+1)*h/radius_segments;

		}

		mode++;
	}
}


// Neighbour Search algorithm function and kernel block
/** Searching for neighbour in particular spatial cell for particular particle
 *  It takes every particles from particular cell and check if it satisfy
 *  a condition that distance between particles is <= closest_distance
 */
void sort(
    float _distanceSquared, 
    __global struct extend_particle * ext_particles, 
    int *last_farthest,
    int neighborParticleId
)
{
    float tmp, tmpId;
    ext_particles->neighbour_list[*last_farthest][0] = neighborParticleId;
    ext_particles->neighbour_list[*last_farthest][1] = _distanceSquared;
    for(int i = *last_farthest; i > 0; --i){
        if(ext_particles->neighbour_list[i - 1][1] > ext_particles->neighbour_list[i][1]){
            tmpId = ext_particles->neighbour_list[i][0];
            tmp = ext_particles->neighbour_list[i][1];
            ext_particles->neighbour_list[i][1] = ext_particles->neighbour_list[i - 1][1];
            ext_particles->neighbour_list[i][0] = ext_particles->neighbour_list[i - 1][0];
            ext_particles->neighbour_list[i - 1][1] = tmp;
            ext_particles->neighbour_list[i - 1][0] = tmpId;
        }
    }
}


//  #define debug_id 1455

void searchForNeighbors(
		int searchCell_,
		__global struct particle * particles,
		__global int * b_grid_cell_id_list,
		float4 position_,
		int myParticleId,
		__global struct extend_particle * ext_particles,
		int *last_farthest,
        float r_thr2
)
{
	if(searchCell_ < 0){
		//printf("\n?????????????? %d ---------- %d ---------- %d\n", searchCell_, grid_offset, myParticleId);
		return ;//last_farthest;
	}
	int baseParticleId = b_grid_cell_id_list[ searchCell_ ];
	if(baseParticleId == -1){
		//printf("\nBAD ASS\n");
		return ;//last_farthest;
	}

	int neighborParticleId = baseParticleId;
	float _distanceSquared;
	// if(myParticleId == debug_id && searchCell_ == 23){
	// 	printf("\n??? %d %d %d ???\n", searchCell_, baseParticleId, particles[neighborParticleId].cell_id);
	// }
	while( particles[neighborParticleId].cell_id == searchCell_){
		if(myParticleId != neighborParticleId)
		{
			float4 d = position_ - particles[neighborParticleId].pos;
			_distanceSquared = d.x * d.x + d.y * d.y + d.z * d.z; // inlined openCL dot(d,d)
			if( _distanceSquared <= r_thr2)
			{
				// if(myParticleId == debug_id && searchCell_ == 23){
				// 	printf("\n??? %d %d ???\n", searchCell_, neighborParticleId);
				// }
                if((*ext_particles).neighbour_list[*last_farthest][1] == -1){
                    sort(_distanceSquared, ext_particles, last_farthest, neighborParticleId);
                    ++(*last_farthest);
                } else if(_distanceSquared < ext_particles->neighbour_list[*last_farthest][1]){
                    sort(_distanceSquared, ext_particles, last_farthest, neighborParticleId);
                }
			}
		}
		++neighborParticleId;
	}
}

/** Searchin for neigbours foe each particles
*/
__kernel void k_neighbour_search(
		__global struct extend_particle * ext_particles,
		__global struct
				particle * particles,
        __global int * b_grid_cell_id_list,
		uint grid_cells_X,
		uint grid_cells_Y,
		uint grid_cells_Z,
		uint grid_cell_count,
        uint grid_offset,
		float h,
		float hashGridCellSize,
		float hashGridCellSizeInv,
		float simulationScale,
		float xmin,
		float ymin,
		float zmin,
		int PARTICLE_COUNT,
		int OFFSET,
		int LIMIT
){
	int id = get_global_id( 0 );
	if(id >= PARTICLE_COUNT || id + OFFSET > LIMIT || id - OFFSET < 0 ){
 		return;
 	}
	float4 position_ = particles[ id ].pos;
	int myCellId = particles[id].cell_id;//& 0xffffff;// truncate to low 16 bits
	int searchCells[8];
	float r_thr2 = h * h;
	float closest_distances[NEIGHBOUR_COUNT];
	int closest_indexes[NEIGHBOUR_COUNT];
	int found_count = 0;
	for(int k=0;k<NEIGHBOUR_COUNT;++k){
		closest_distances[k] = r_thr2;
		closest_indexes[k] = -1;
	}
	searchCells[0] = myCellId ;//- grid_offset;
	// p is the current particle position within the bounds of the hash grid
	float4 p;
	float4 p0 = (float4)( xmin, ymin, zmin, 0.0f );
	p = position_ - p0;

	// cf is the min,min,min corner of the current cell
	int4 cellFactors_ = cellFactors( position_, xmin, ymin, zmin, hashGridCellSizeInv );
	float4 cf;
	cf.x = cellFactors_.x * hashGridCellSize;
	cf.y = cellFactors_.y * hashGridCellSize;
	cf.z = cellFactors_.z * hashGridCellSize;

	// lo.A is true if the current position is in the low half of the cell for dimension A
	//float4 ttt = ( p - cf );
	int4 lo;
	
	lo = isless(p, h + cf);
	int p_id = id;//particles[id].particle_id;
	int4 delta;
    int4 one = (int4)( 1, 1, 1, 1 );
	delta = one + 2 * lo;
	int last_farthest = 0;

	searchCells[1] = searchCell( myCellId, delta.x, 0,       0,       grid_cells_X, grid_cells_Y, grid_cells_Z, grid_cell_count ) ;//- grid_offset;
	searchCells[2] = searchCell( myCellId, 0,       delta.y, 0,       grid_cells_X, grid_cells_Y, grid_cells_Z, grid_cell_count ) ;//- grid_offset;
	searchCells[3] = searchCell( myCellId, 0,       0,       delta.z, grid_cells_X, grid_cells_Y, grid_cells_Z, grid_cell_count ) ;//- grid_offset;
	searchCells[4] = searchCell( myCellId, delta.x, delta.y, 0,       grid_cells_X, grid_cells_Y, grid_cells_Z, grid_cell_count ) ;//- grid_offset;
	searchCells[5] = searchCell( myCellId, delta.x, 0,       delta.z, grid_cells_X, grid_cells_Y, grid_cells_Z, grid_cell_count ) ;//- grid_offset;
	searchCells[6] = searchCell( myCellId, 0,       delta.y, delta.z, grid_cells_X, grid_cells_Y, grid_cells_Z, grid_cell_count ) ;//- grid_offset;
	searchCells[7] = searchCell( myCellId, delta.x, delta.y, delta.z, grid_cells_X, grid_cells_Y, grid_cells_Z, grid_cell_count ) ;//- grid_offset;

	//  if(p_id == debug_id) {
	//  	printf("\n%d, %d, %d, %d, %d, %d, %d, %d\n", searchCells[0], searchCells[1], searchCells[2], searchCells[3],
	//  	       searchCells[4], searchCells[5], searchCells[6], searchCells[7]);
	//  }
	// Search neighbour particles in every cells from searchCells list
	searchForNeighbors( searchCells[0], particles, b_grid_cell_id_list, position_,
	                                      id, &ext_particles[id], &last_farthest, r_thr2 );

	searchForNeighbors( searchCells[1], particles, b_grid_cell_id_list, position_,
	                                      id, &ext_particles[id], &last_farthest, r_thr2 );

	searchForNeighbors( searchCells[2], particles, b_grid_cell_id_list, position_,
	                                      id, &ext_particles[id], &last_farthest, r_thr2 );


	searchForNeighbors( searchCells[3], particles, b_grid_cell_id_list, position_,
	                                      id, &ext_particles[id], &last_farthest, r_thr2 );


	searchForNeighbors( searchCells[4], particles, b_grid_cell_id_list, position_,
	                                      id, &ext_particles[id], &last_farthest, r_thr2 );

	searchForNeighbors( searchCells[5], particles, b_grid_cell_id_list, position_,
	                                      id, &ext_particles[id], &last_farthest, r_thr2 );

	searchForNeighbors( searchCells[6], particles, b_grid_cell_id_list, position_,
	                                      id, &ext_particles[id], &last_farthest, r_thr2 );

	searchForNeighbors( searchCells[7], particles, b_grid_cell_id_list, position_,
	                                      id, &ext_particles[id], &last_farthest, r_thr2 );

    // storing all found neighbors into neighborMap buffer
	for(int j=0; j<NEIGHBOUR_COUNT; ++j){
        if((int)ext_particles[id].neighbour_list[j][0] == NO_PARTICLE_ID)
            break;
		ext_particles[id].neighbour_list[j][1] = SQRT( ext_particles[id].neighbour_list[j][1] ) * simulationScale;
		if(ext_particles[id].neighbour_list[j][1] == 0.0f) {
	 	 	printf("\n======\tP_ID: %d, CELL:_ID: %d,\n", 
		 	 (int)ext_particles[id].neighbour_list[j][0], 
		 	 particles[(int)ext_particles[id].neighbour_list[j][0]].cell_id
		 	 );
	 	}
	}
}