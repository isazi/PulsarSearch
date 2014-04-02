// Copyright 2013 Alessio Sclocco <a.sclocco@vu.nl>
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <string>
using std::string;
#include <map>
using std::map;
#include <vector>
using std::vector;


#ifndef READ_CONFIGURATION_HPP
#define READ_CONFIGURATION_HPP

void readPadding(map< string, unsigned int > & padding, const string paddingFilename);
void readVectorWidth(map< string, unsigned int > & vectorWidth, const string vectorFilename);
void readDedispersion(map< string, map< unsigned int, vector< unsigned int > > > & dedispersionParameters, const string dedispersionFilename);
void readTranspose(map< string, map< unsigned int, unsigned int > > & transposeParameters, const string transposeFilename);
void readFolding(std::map< std::string, std::map< unsigned int, std::map< unsigned int, std::vector< unsigned int > > > > & foldingParameters, const string foldingFilename);
void readSNR(std::map< std::string, std::map< unsigned int, std::map< unsigned int, std::vector< unsigned int > > > > & snrParameters, const string snrFilename);


#endif // READ_CONFIGURATION_HPP
