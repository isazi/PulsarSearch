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
#include <map>
#include <vector>


#ifndef READ_CONFIGURATION_HPP
#define READ_CONFIGURATION_HPP

void readPadding(std::map< std::string, unsigned int > & padding, const std::string & paddingFilename);
void readVectorWidth(std::map< std::string, unsigned int > & vectorWidth, const std::string & vectorFilename);
void readDedispersion(std::map< std::string, std::map< unsigned int, std::vector< unsigned int > > > & dedispersionParameters, const std::string & dedispersionFilename);
void readTranspose(std::map< std::string, std::map< unsigned int, unsigned int > > & transposeParameters, const std::string & transposeFilename);
void readFolding(std::std::map< std::std::string, std::std::map< unsigned int, std::std::map< unsigned int, std::vector< unsigned int > > > > & foldingParameters, const std::string & foldingFilename);
void readSNR(std::std::map< std::std::string, std::std::map< unsigned int, std::std::map< unsigned int, std::vector< unsigned int > > > > & snrParameters, const std::string & snrFilename);


#endif // READ_CONFIGURATION_HPP

