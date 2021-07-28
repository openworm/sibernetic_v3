#include <utility>

/*******************************************************************************
 * The MIT License (MIT)
 *
 * Copyright (c) 2011, 2017 OpenWorm.
 * http://openworm.org
 *
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the MIT License
 * which accompanies this distribution, and is available at
 * http://opensource.org/licenses/MIT
 *
 * Contributors:
 *     	OpenWorm - http://openworm.org/people.html
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
 * USE OR OTHER DEALINGS IN THE SOFTWARE.
 *******************************************************************************/
#ifndef OCLSOLVER_HPP
#define OCLSOLVER_HPP

#if defined(_WIN32) || defined(_WIN64)
#pragma comment(lib, "opencl.lib") // opencl.lib
#endif

#if defined(__APPLE__) || defined(__MACOSX)

#include "OpenCL/cl.hpp"

#else

#include <CL/cl.hpp>

#endif

#include "isolver.h"
#include "ocl_const.h"
#include "ocl_struct.h"
#include "particle.h"
#include "sph_model.hpp"
#include "util/error.h"
#include "util/json.hpp"
#include "device.h"
#include <fstream>
#include <iostream>

namespace sibernetic {
	namespace solver {
		using std::cout;
		using std::endl;
		using std::shared_ptr;
		using sibernetic::model::particle;
		using sibernetic::model::partition;
		using sibernetic::model::sph_model;
		using sibernetic::ocl_error;
		enum LOGGING_MODE {
			FULL,
			NO
		};
// OCL constans block
#define QUEUE_EACH_KERNEL 1

		const int LOCAL_NDRANGE_SIZE = 256;

		template<class T = float>
		class ocl_solver : public i_solver {
			typedef shared_ptr<sph_model<T>> model_ptr;

		public:
			ocl_solver(
					model_ptr &m,
					shared_ptr<device> d,
					size_t idx,
					LOGGING_MODE log_mode = LOGGING_MODE::NO):
				model(m),
				dev(std::move(d)),
				device_index(idx),
				log_mode(log_mode)
			{
				try {
					this->initialize_ocl();
				} catch (ocl_error &ex) {
					throw;
				}
			}
            std::shared_ptr<device> get_device() override {
                return this->dev;
			}
			// TODO rename method!!!
			void init_model(partition *p) override {
				this->p = p;
				prev_part_size = p->total_size();
				init_buffers();
				init_kernels();
			}

			~ocl_solver() override = default;

			void _debug_(){
				std::vector<extend_particle> neighbour_map(p->size());
				copy_buffer_from_device(&(neighbour_map[0]), b_ext_particles, p->size() * sizeof(extend_particle), 0);
				copy_buffer_from_device(&(model->get_particles()[0]), b_particles, p->size() * sizeof(particle<T>), 0);
				std::string big_s = "[";
				for(auto p: neighbour_map){

					big_s += "{\"particle\": ";
					big_s += model->get_particle(p.p_id).jsonify();
					big_s += ",";
//					big_s += "\"particle_id\": ";
//					big_s += std::to_string(p.p_id);
//					big_s += ",";
					big_s += "\"n_list\":[";
					for(int i=0;i<NEIGHBOUR_COUNT;++i) {
						big_s += "{";
						big_s += "\"n_particle_id\": ";
						big_s += std::to_string(p.neighbour_list[i][0]);
						big_s += ",";
						big_s += "\"distance\": ";
						big_s += std::to_string(p.neighbour_list[i][1]);
						big_s += "}";
						if(i != NEIGHBOUR_COUNT - 1)
							big_s += ",";
					}
					big_s += "]";
					big_s += "}";
					//break;
					if(p.p_id != neighbour_map.back().p_id)
						big_s += ",";
				}
				big_s += "]";
				std::ofstream debug_file("debug");
				debug_file << big_s;
				debug_file.close();
			}

			void neighbour_search() override {
				run_init_ext_particles();
				// run_hash_particles();
				// sync();
				// run_clear_grid_hash();
				// run_fill_particle_cell_hash();
				// run_neighbour_search();
			}

			void physic() override {
				// int iter = 0;
				// run_compute_density();
				// run_compute_forces_init_pressure();
				// while(iter < sibernetic::model::PCI_ITER_COUNT) {
				// 	run_predict_positions();
				// 	run_predict_density();
				// 	run_correct_pressure();
				// 	run_compute_pressure_force_acceleration();
				// 	++iter;
				// }
				// run_integrate();
			}

			void sync() override {
				is_synchronizing = true;
				std::cout << "start " <<  p->start << " end " <<  p->end << " size = " << p->size() << " offset = " << p->offset() * sizeof(particle<T>) << " len of buff host " << model->get_particles().size() << std::endl;
				copy_buffer_from_device(
						&(model->get_particles()[p->start]),
						b_particles,
						p->size() * sizeof(particle<T>),
						p->offset() * sizeof(particle<T>));
				if(model->set_ready()){
					model->sync();
				} else {
					while(is_synchronizing);
				}

				init_buffers();
				//copy_buffer_to_device((void *) &(model->get_particles()[p->ghost_start]), b_particles, 0, p->total_size() * sizeof(particle<T>));
			}

		void run(int iter_lim) override {
			int i = 0;
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wmissing-noreturn"
			while(true) {
				if(iter_lim != -1 && i == iter_lim) {
					//_debug_();
					break;
				}
				neighbour_search();
				physic();
				++i;
			}
#pragma clang diagnostic pop
		}
		void unfreeze() override {
			is_synchronizing = false;
		}
		private:
			std::atomic<bool> is_synchronizing;
			model_ptr model;
			int prev_part_size;
			int prev_start;
			int prev_end;
			size_t device_index;
			partition* p;
			shared_ptr<device> dev;
			std::string msg = dev->name + '\n';
			const std::string cl_program_file = "cl_code//sph_cl_code.cl";
			LOGGING_MODE log_mode;
			cl::Kernel k_init_ext_particles;
			cl::Kernel k_hash_particles;
			cl::Kernel k_clear_grid_hash;
			cl::Kernel k_fill_particle_cell_hash;
			cl::Kernel k_neighbour_search;

			cl::Kernel k_compute_density;
			cl::Kernel k_compute_forces_init_pressure;
			cl::Kernel k_predict_positions;
			cl::Kernel ker_predict_density;
			cl::Kernel ker_correct_pressure;
			cl::Kernel ker_compute_pressure_force_acceleration;
			cl::Kernel k_integrate;

			cl::Buffer b_particles;
			cl::Buffer b_ext_particles;
			cl::Buffer b_grid_cell_id_list;
			cl::CommandQueue queue;
			cl::Program program;

			void init_buffers() {
				auto p_buff_size = create_ocl_buffer("particles", b_particles, CL_MEM_READ_WRITE,
				                  p->total_size() * sizeof(particle<T>));
				auto ext_p_buff_size = create_ocl_buffer("ext_particles", b_ext_particles, CL_MEM_READ_WRITE,
				                  p->total_size() * sizeof(extend_particle));
				auto b_grid_buff_size = create_ocl_buffer("b_grid_cell_id_list", b_grid_cell_id_list, CL_MEM_READ_WRITE,
				                  model->get_total_cell_num() * sizeof(int));
				copy_buffer_to_device((void *) &(model->get_particles()[p->ghost_start]),
				                      b_particles, 0, p->total_size() * sizeof(particle<T>));
				auto total_device_alloc_memory = (p_buff_size + ext_p_buff_size + b_grid_buff_size) / (1024 * 1024);
				auto ram_memory = (p->total_size() * sizeof(particle<T>) + p->total_size() * sizeof(extend_particle) + model->get_total_cell_num() * sizeof(int))/(1024 * 1024);
				//std::cout << "Size of object on RAM " << ram_memory << "MB.Size of object on DEVICE RAM " << total_device_alloc_memory  << "MB."<< std::endl;
			}

			void init_kernels() {
				create_ocl_kernel("k_init_ext_particles", k_init_ext_particles);
				create_ocl_kernel("k_hash_particles", k_hash_particles);
				create_ocl_kernel("k_clear_grid_hash", k_clear_grid_hash);
				create_ocl_kernel("k_fill_particle_cell_hash", k_fill_particle_cell_hash);
				create_ocl_kernel("k_neighbour_search", k_neighbour_search);

				create_ocl_kernel("k_compute_density", k_compute_density);
				create_ocl_kernel("k_compute_forces_init_pressure", k_compute_forces_init_pressure);
				create_ocl_kernel("k_predict_positions", k_predict_positions);
				create_ocl_kernel("ker_predict_density", ker_predict_density);
				create_ocl_kernel("ker_correct_pressure", ker_correct_pressure);
				create_ocl_kernel("ker_compute_pressure_force_acceleration", ker_compute_pressure_force_acceleration);
				create_ocl_kernel("k_integrate", k_integrate);
			}

			void init_ext_particles() override {}

			void initialize_ocl() {
				int err;
				queue = cl::CommandQueue(dev->context, dev->dev, 0, &err);
				if (err != CL_SUCCESS) {
					throw ocl_error(msg + "Failed to create command queue");
				}
				std::ifstream file(cl_program_file);
				if (!file.is_open()) {
					throw ocl_error(msg + "Could not open file with OpenCL program check "
					                      "input arguments oclsourcepath: " +
					                cl_program_file);
				}
				std::string programSource(std::istreambuf_iterator<char>(file),
				                          (std::istreambuf_iterator<char>()));
				// TODO fix this to param
				if (0) {
					programSource =
							"#define _DOUBLE_PRECISSION\n" +
							programSource; // not now it needs double extension check on device
				}
				cl::Program::Sources source(
						1, std::make_pair(programSource.c_str(), programSource.length() + 1));
				program = cl::Program(dev->context, source);
#if defined(__APPLE__)
				err = program.build("-g -cl-opt-disable -I .");
#else
#if INTEL_OPENCL_DEBUG
				err = program.build(OPENCL_DEBUG_PROGRAM_PATH + "-g -cl-opt-disable -I .");
#else
				err = program.build("-I .");
#endif
#endif
				if (err != CL_SUCCESS) {
					std::string compilationErrors;
					compilationErrors = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(dev->dev);
					msg += make_msg(msg, "Compilation failed: ", compilationErrors);
					throw ocl_error(msg);
				}
				std::cout
						<< msg
						<< "OPENCL program was successfully build. Program file oclsourcepath: "
						<< cl_program_file << std::endl;
			}

            size_t create_ocl_buffer(const char *name, cl::Buffer &b,
			                       const cl_mem_flags flags, const size_t size) {
				int err;
				b = cl::Buffer(dev->context, flags, size, nullptr, &err);
				if (err != CL_SUCCESS) {
					std::string error_m =
							make_msg("Buffer creation failed: ", name, " Error code is ", err);
					throw ocl_error(error_m);
				}
				size_t ocl_obj_size;
				b.getInfo(CL_MEM_SIZE, &ocl_obj_size);
				//std::cout << "Size of buffer " << name << " on RAM " << size / (1024 * 1024) << "MB. Size of object on DEVICE RAM " << ocl_obj_size / (1024 * 1024) << " MB"<< std::endl;
				return ocl_obj_size;
			}
			void create_ocl_kernel(const char *name, cl::Kernel &k) {
				int err;
				k = cl::Kernel(program, name, &err);
				if (err != CL_SUCCESS) {
					std::string error_m =
							make_msg("Kernel creation failed: ", name, " Error code is ", err);
					throw ocl_error(error_m);
				}
			}

			void copy_buffer_to_device(
					const void *host_b,
					cl::Buffer &ocl_b,
			        const size_t offset,
			        const size_t size) {
				// Actually we should check  size and type
				int err = queue.enqueueWriteBuffer(ocl_b, CL_TRUE, offset, size, host_b);

				if (err != CL_SUCCESS) {
					std::string error_m =
							make_msg("Copy buffer to device is failed error code is ", err);
					throw ocl_error(error_m);
				}
				queue.finish();
			}

			void copy_buffer_from_device(void *host_b, const cl::Buffer &ocl_b,
			                             const size_t size, size_t offset) {
				// Actualy we should check  size and type
				//std::cout << "size = " << size << " offset = " << offset << std::endl;
				int err = queue.enqueueReadBuffer(ocl_b, CL_TRUE, offset, size, host_b);

				if (err != CL_SUCCESS) {
					std::string error_m =
							make_msg("Copy buffer from device is failed error code is ", err);
					throw ocl_error(error_m);
				}
				queue.finish();
			}

			int run_init_ext_particles() {
				if(log_mode == LOGGING_MODE::FULL)
					std::cout << "run init_ext_particles --> " << dev->name << std::endl;
				this->kernel_runner(
						this->k_init_ext_particles,
						p->total_size(),
						0,
						this->b_ext_particles,
						p->total_size()
				);
			}

			int run_hash_particles() {
				if(log_mode == LOGGING_MODE::FULL)
					std::cout << "run hash_particles --> " << dev->name << std::endl;
				this->kernel_runner(
						this->k_hash_particles,
						p->total_size(),
						0,
						this->b_particles,
						model->get_cell_num_x(),
						model->get_cell_num_y(),
						model->get_cell_num_z(),
						sibernetic::model::GRID_CELL_SIZE_INV,
						p->total_size(),
						p->offset(),
						p->limit()
				);
			}

			int run_clear_grid_hash() {
				if(log_mode == LOGGING_MODE::FULL)
					std::cout << "run clear_grid_hash --> " << dev->name << std::endl;

				this->kernel_runner(
						this->k_clear_grid_hash,
						model->get_total_cell_num(),//p.total_cell_count(),
						0,
						this->b_grid_cell_id_list,
						model->get_total_cell_num()
				);
			}

			int run_fill_particle_cell_hash() {
				if(log_mode == LOGGING_MODE::FULL)
					std::cout << "run fill_particle_cell_hash --> " << dev->name << std::endl;

				this->kernel_runner(
						this->k_fill_particle_cell_hash,
						p->total_size(),
						0,
						this->b_grid_cell_id_list,
						this->b_particles,
						p->start_ghost_cell_id,
						p->total_size()
				);
			}

			int run_neighbour_search() {
				if(log_mode == LOGGING_MODE::FULL)
					std::cout << "run neighbour_search --> " << dev->name << std::endl;

				this->kernel_runner(
						this->k_neighbour_search,
						p->total_size(),
						0,
						this->b_ext_particles,
						this->b_particles,
						this->b_grid_cell_id_list,
						model->get_cell_num_x(),
						model->get_cell_num_y(),
						model->get_cell_num_z(),
						model->get_total_cell_num(),
						p->start_ghost_cell_id,
						sibernetic::model::H,
						sibernetic::model::GRID_CELL_SIZE,
						sibernetic::model::GRID_CELL_SIZE_INV,
						model->get_config()["simulation_scale"],
						model->get_config()["x_min"],
						model->get_config()["y_min"],
						model->get_config()["z_min"],
						p->total_size(),
						p->offset(),
						p->limit()
				);
			}

			int run_compute_density() {
				if(log_mode == LOGGING_MODE::FULL)
					std::cout << "run compute_density --> " << dev->name << std::endl;

				this->kernel_runner(
						this->k_compute_density,
						p->total_size(),
						0,
						this->b_ext_particles,
						this->b_particles,
						model->get_config()["mass_mult_wpoly6_coefficient"],
						model->get_config()["h_scaled_2"],
						p->total_size(),
						p->offset(),
						p->limit()
				);
			}

			int run_compute_forces_init_pressure() {
				if(log_mode == LOGGING_MODE::FULL)
					std::cout << "run compute_forces_init_pressure --> " << dev->name << std::endl;

				this->kernel_runner(
						this->k_compute_forces_init_pressure,
						p->total_size(),
						0,
						this->b_ext_particles,
						this->b_particles,
						model->get_config()["surf_tens_coeff"],
						model->get_config()["mass_mult_divgrad_viscosity_coefficient"],
						model->get_config()["h_scaled"],
						model->get_config()["gravity_x"],
						model->get_config()["gravity_y"],
						model->get_config()["gravity_z"],
						p->total_size(),
						p->offset(),
						p->limit()
				);
			}

			void run_predict_positions(){
				if(log_mode == LOGGING_MODE::FULL)
					std::cout << "run predict_positions --> " << dev->name << std::endl;

				this->kernel_runner(
						this->k_predict_positions,
						p->total_size(),
						0,
						this->b_ext_particles,
						this->b_particles,
						model->get_config()["simulation_scale_inv"],
						model->get_config()["time_step"],
						sibernetic::model::R_0,
						p->total_size(),
						p->offset(),
						p->limit()
				);
			}

			void run_predict_density(){
				if(log_mode == LOGGING_MODE::FULL)
					std::cout << "run predict_density --> " << dev->name << std::endl;

				this->kernel_runner(
						this->ker_predict_density,
						p->total_size(),
						0,
						this->b_ext_particles,
						this->b_particles,
						model->get_config()["mass_mult_wpoly6_coefficient"],
						sibernetic::model::H,
						model->get_config()["simulation_scale"],
						p->total_size(),
						p->offset(),
						p->limit()
				);
			}
			void run_correct_pressure(){
				if(log_mode == LOGGING_MODE::FULL)
					std::cout << "run run_correct_pressure --> " << dev->name << std::endl;

				this->kernel_runner(
						this->ker_correct_pressure,
						p->total_size(),
						0,
						this->b_particles,
						sibernetic::model::DENSITY_WATER,
						model->get_config()["delta"],
						p->total_size(),
						p->offset(),
						p->limit()
				);
			}
			void run_compute_pressure_force_acceleration(){
				if(log_mode == LOGGING_MODE::FULL)
					std::cout << "run run_compute_pressure_force_acceleration --> " << dev->name << std::endl;

				this->kernel_runner(
						this->ker_compute_pressure_force_acceleration,
						p->total_size(),
						0,
						this->b_ext_particles,
						this->b_particles,
						model->get_config()["mass_mult_grad_wspiky_coefficient"],
						model->get_config()["h_scaled"],
						model->get_config()["simulation_scale"],
						model->get_config()["delta"],
						sibernetic::model::DENSITY_WATER,
						p->total_size(),
						p->offset(),
						p->limit()
				);
			}
			void run_integrate(){
				if(log_mode == LOGGING_MODE::FULL)
					std::cout << "run run_integrate --> " << dev->name << std::endl;

				this->kernel_runner(
						this->k_integrate,
						p->total_size(),
						0,
						this->b_ext_particles,
						this->b_particles,
						model->get_config()["simulation_scale_inv"],
						model->get_config()["time_step"],
						sibernetic::model::R_0,
						1,
						2,
						p->total_size(),
						p->offset(),
						p->limit()
				);
			}

			template<typename U, typename... Args>
			int kernel_runner(
					cl::Kernel &ker,
					unsigned int dim,
					unsigned int arg_pos,
					U &arg,
					Args... args) {
				ker.setArg(arg_pos, arg);
				return kernel_runner(ker, dim, ++arg_pos, args...);
			}

			template<typename U>
			int kernel_runner(cl::Kernel &ker, unsigned int dim, unsigned int arg_pos, U &arg) {
			    auto lk = this->dev->global_work_group_size;
				ker.setArg(arg_pos, arg);
				auto dim_round_up = (((dim - 1) / lk) + 1) * lk;
				auto err = queue.enqueueNDRangeKernel(
						ker, cl::NullRange, cl::NDRange(dim_round_up),
#if defined(__APPLE__)
						cl::NullRange, nullptr, nullptr);
#else
						cl::NDRange(this->dev->global_work_group_size), nullptr, nullptr);
#endif
#if QUEUE_EACH_KERNEL
				queue.finish();
#endif
				if (err != CL_SUCCESS) {
					std::string k_name;
					ker.getInfo(CL_KERNEL_FUNCTION_NAME, &k_name);
					throw ocl_error(make_msg("Kernel finish its work with errror", err, "kernel name",k_name));
				}
				return err;
			}

		};
	} // namespace solver
} // namespace sibernetic
#endif // OCLSOLVER_HPP
