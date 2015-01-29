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

#include <Dedispersion.hpp>
#include <Transpose.hpp>

#ifndef READ_CONFIGURATION_HPP
#define READ_CONFIGURATION_HPP

void readPadding(std::map< std::string, unsigned int > & padding, const std::string & paddingFilename);
void readVectorWidth(std::map< std::string, unsigned int > & vectorWidth, const std::string & vectorFilename);
void readDedispersion(std::map< std::string, std::map< unsigned int, PulsarSearch::DedispersionConf > > & dedispersionParameters, const std::string & dedispersionFilename);
void readTranspose(std::map< std::string, std::map< unsigned int, isa::OpenCL::transposeConf > > & transposeParameters, const std::string & transposeFilename);
void readSNRD(std::map< std::string, std::map< unsigned int, PulsarSearch::snrDedispersedConf > > & snrParameters, const std::string & snrFilename);
void readFolding(std::map< std::string, std::map< unsigned int, std::map< unsigned int, PulsarSearch::FoldingConf > > > & foldingParameters, const std::string & foldingFilename);
void readSNRF(std::map< std::string, std::map< unsigned int, std::map< unsigned int, PulsarSearch::snrFoldedConf > > > & snrParameters, const std::string & snrFilename);

#endif // READ_CONFIGURATION_HPP

