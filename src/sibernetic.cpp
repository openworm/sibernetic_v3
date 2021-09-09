/*******************************************************************************
 * The MIT License (MIT)
 *
 * Copyright (c) 2011, 2013 OpenWorm.
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
#include "solver_container.hpp"
#include "isolver.h"
#include "util/arg_parser.h"
#include "util/custom_reader.hpp"
#include <iostream>
#include <thread>
#include <time.h>

using sibernetic::model::sph_model;
using sibernetic::solver::solver_container;
using sibernetic::solver::EXECUTION_MODE;


int main(int argc, char **argv) {
  arg_parser prsr(argc, argv);
  if (prsr.check_arg("-h") || prsr.check_arg("--help") ||
      prsr.check_arg("-?") || prsr.check_arg("-help")) {
    return arg_parser::show_usage();
  }
  std::string model_name;
  auto mode = EXECUTION_MODE::ONE;
  bool parallel_sort = false;
  size_t dev_count = 0;
  if (prsr.check_arg("-f")) {
    model_name = prsr.get_arg_value("-f");
  } else {
    model_name = "config/data/1.txt";
  }
  if (prsr.check_arg("--multi_dev")) {
  	auto v = prsr.get_arg_value("--multi_dev");
  	if(v != "ALL"){
		dev_count = std::stol(v);
  	}
  	mode = EXECUTION_MODE::ALL;
  }
  if (prsr.check_arg("--p_sort")) {
    parallel_sort = true;
  }
  int iter_lim = -1;
  if (prsr.check_arg("-iter_lim")) {
    iter_lim = std::atoi(prsr.get_arg_value("-iter_lim").c_str());
  }

  float time_lim = -1;
  if (prsr.check_arg("-time_lim")) {
    time_lim = std::atof(prsr.get_arg_value("-time_lim").c_str());
  }
  try {
    abstract_reader<float> * reader = nullptr;
    if (prsr.check_arg("--json")) {
        reader = new json_reader<float>();
    } else {
        reader = new custom_reader<float>();
    }

    std::shared_ptr<sph_model<float>> model(new sph_model<float>(model_name, reader));
    solver_container<float> &s_con = solver_container<float>::instance(
            model,
            mode,
            dev_count,
            sibernetic::solver::OCL,
            parallel_sort
    );
    //model->sort();
    timespec t1, t2;
    clock_gettime(CLOCK_MONOTONIC_RAW, &t1);
    s_con.run(iter_lim, time_lim);
    clock_gettime(CLOCK_MONOTONIC_RAW, &t2);
    auto delta = ((t2.tv_sec - t1.tv_sec) * 1000000 + (t2.tv_nsec - t1.tv_nsec) / 1000 ) / 1000000.0 ;
    std::cout << "Simulation has run " << delta << " s" << std::endl;
  } catch (sibernetic::parser_error &e) {
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (sibernetic::ocl_error &e) {
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
