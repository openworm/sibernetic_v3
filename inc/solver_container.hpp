#ifndef SOLVER_CONTAINER_HPP
#define SOLVER_CONTAINER_HPP

#include "isolver.h"
#include "ocl_const.h"
#include "ocl_solver.hpp"
#include "ocl_sort_solver.hpp"
#include "ocl_radix_sort.hpp"
#include "sph_model.hpp"
#include "util/ocl_helper.h"
#include "util/error.h"
#include <string>
#include <vector>
#include <thread>
#include <cmath>
#include <exception>

namespace sibernetic {
	namespace solver {
		using model::sph_model;
		using std::shared_ptr;
		using sibernetic::solver::ocl_solver;
        using sibernetic::solver::ocl_sort_solver;
        using sibernetic::solver::ocl_radix_sort_solver;
		enum EXECUTION_MODE {
			ONE,
			ALL
		};

		template<class T = float>
		class solver_container {
			typedef shared_ptr<sph_model<T>> model_ptr;

		public:
			solver_container(const solver_container &) = delete;

			solver_container &operator=(const solver_container &) = delete;

			/** Maer's singleton
			 */
			static solver_container &instance(
			        model_ptr &model,
					EXECUTION_MODE mode = EXECUTION_MODE::ONE,
					size_t dev_count = 0,
			        SOLVER_TYPE s_t = OCL,
			        bool p_sort = false
			    ) {
				static solver_container s(model, mode, s_t, dev_count, p_sort);
				model->set_solver_container(s.solvers());
				return s;
			}

			void run(int iter_lim, float time_lim) {
				std::vector<std::thread> t_pool;
				std::for_each(
					_solvers.begin(), 
					_solvers.end(), 
					[&, this](std::shared_ptr<i_solver> &s) {
						t_pool.emplace_back(std::thread(solver_container::run_solver, std::ref(s), iter_lim, time_lim));
					}
				);
				std::for_each(t_pool.begin(), t_pool.end(), [](std::thread &t) { t.join(); });
				if(finished_with_error){
					std::cout << "when runnint ocl solver some problem was occured store fail config into debug.out file" << std::endl;
					for(auto log: logs){
						std::cout << *log << std::endl;
					}
					_debug_();
					return;
				}
			}

			static void run_solver(std::shared_ptr<i_solver> &s, int iter_lim, float time_lim) {
				s->run(iter_lim, time_lim);
			}
			std::vector<std::shared_ptr<i_solver>>* solvers(){
				return &_solvers;
			}
		private:
			explicit solver_container(model_ptr &model, EXECUTION_MODE mode = EXECUTION_MODE::ONE,
			                          SOLVER_TYPE s_type = OCL, size_t dev_count = 0, bool p_sort = false) {
				try {
					p_q dev_q = get_dev_queue();
					size_t device_index = 0;
					while (!dev_q.empty()) {
						try {
							std::string * log = new(std::string);
							logs.push_back(log);
							std::shared_ptr<ocl_solver<T>> solver(new ocl_solver<T>(model, dev_q.top(), device_index, &finished_with_error, log));
							_solvers.push_back(solver);
							std::cout << "************* DEVICE For Phys *************" << std::endl;
							dev_q.top()->show_info();
							std::cout << "**********************************" << std::endl;
							++device_index;
                            dev_q.pop();
							if(mode == EXECUTION_MODE::ONE){
								break;
							} else if(mode == EXECUTION_MODE::ALL){
								if(dev_count > 0 && dev_count == device_index){
								    break;
								}
							}
						} catch (ocl_error &ex) {
							std::cout << ex.what() << std::endl;
						}
					}

					if(!dev_q.empty() && p_sort) {
					    std::cout << "************* DEVICE 4 SORT *************" << std::endl;
                        dev_q.top()->show_info();
                        std::cout << "*******************************************" << std::endl;

					    std::shared_ptr<ocl_radix_sort_solver<T>> sort_solver(new ocl_radix_sort_solver<T>(model, dev_q.top(), device_index+1));
					    sort_solver->init_model();
					    model->set_sort_solver(sort_solver);
					}

					if (_solvers.size()) {
			                        init_weights();
                        			model->set_balance_vector(this->weights);
						model->make_partition(_solvers.size()); // TODO to think about is in future we
                        //model->make_partition(_solvers.size(), std::vector<float>{0.5, 0.4, 0.1});
						// can't init one or more
						// devices
						// obvious we should reinit partitions case ...
						int i=0;
						for (auto s : _solvers) {
						    s->get_device()->balance_coeff = weights[i++];
							s->init_model(&(model->get_next_partition()));
						}
					} else {
						throw ocl_error("No OpenCL devices were initialized.");
					}
					this->md = model;
				} catch (sibernetic::ocl_error &err) {
					throw;
				}
			}

            void init_weights(){
			    size_t total_compute_power = 0;
			    std::for_each(_solvers.begin(), _solvers.end(),
			            [&total_compute_power](std::shared_ptr<i_solver> s) -> void {
			                total_compute_power += s->get_device()->max_thread_count;
			        }
			    );
                float rest = 1.f;
			    for(int i=0; i < _solvers.size(); ++i){
                    float cur;
			        if(i == _solvers.size() - 1) {
			            cur = ceilf(rest * 100) / 100;
			        }else {
			            auto s = _solvers[i];
			            cur = 1.f / static_cast<float>(_solvers.size()); //static_cast<float>(s->get_device()->max_thread_count) / static_cast<float>(total_compute_power);
			            cur = (cur * 100 + 0.5f) / 100.f;
			            rest -= cur;
			        }
			       weights.push_back(cur);
				//weights.push_back(0.3f);
			    }
			}

			~solver_container() = default;
            std::vector<float> weights;
			model_ptr md;
			std::vector<std::shared_ptr<i_solver>> _solvers;
            std::shared_ptr<i_sort_solver> sort_solver;
			bool finished_with_error;
			std::vector<std::string*> logs;
			void _debug_(){
				std::vector<extend_particle> neighbour_map(md->size());
				
				std::stringstream big_s;
				big_s << "bbox\n";
				big_s << md->get_config()["x_max"];
				big_s << "\n";
				big_s << md->get_config()["x_min"];
				big_s << "\n";
				big_s << md->get_config()["y_max"];
				big_s << "\n";
				big_s << md->get_config()["y_min"];
				big_s << "\n";
				big_s << md->get_config()["z_max"];
				big_s << "\n";
				big_s << md->get_config()["z_min"];
				big_s << "\n";
				big_s << "position\n";
				for(auto p: md->get_particles()) {
					big_s << p.pos[0];
					big_s << "\t";
					big_s << p.pos[1];
					big_s << "\t";
					big_s << p.pos[2];
					big_s << "\t";
					big_s << p.pos[3];
					big_s << "\t";
					big_s << p.vel[0];
					big_s << "\t";
					big_s << p.vel[1];
					big_s << "\t";
					big_s << p.vel[2];
					big_s << "\t";
					big_s << p.vel[3];
					big_s << "\n";
				}
				
				std::ofstream debug_file("debug.out");
				debug_file << big_s.str();
				debug_file.close();
			}
		};
	} // namespace solver
} // namespace sibernetic

#endif // SOLVER_CONTAINER_HPP
