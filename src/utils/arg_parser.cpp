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
#include "util/arg_parser.h"
#include <algorithm>

arg_parser::arg_parser(int argc, char **argv) {
  for (int i = 1; i < argc; ++i) {
  	auto s = std::string(argv[i]);
    arguments.insert({s.substr(0,s.find('=')), s});//.push_back(s.substr(s.find('=')));
  }
}

bool arg_parser::check_arg(const std::string &arg) const {
  return arguments.find(arg.substr(0,arg.find('='))) != arguments.end();
}

const std::string arg_parser::get_arg_value(const std::string &arg) const {
	auto itr = arguments.find(arg);
	if (itr != arguments.end()) {
		auto val = arguments.at(arg);
		if (val.find('=') != std::string::npos) {
			return val.substr(val.find('=') + 1);
		}
	}

	return std::string("ALL");
}
