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
#include <fstream>
#include <map>
#include <vector>
#include <cctype>

#include <readConfiguration.hpp>
#include <utils.hpp>


void readPadding(map< std::string, unsigned int > & padding, const std::string & paddingFilename) {
	std::string temp;
	std::ifstream paddingFile(paddingFilename);

	while ( ! paddingFile.eof() ) {
		unsigned int middle = 0;

		std::getline(paddingFile, temp);
		if ( ! std::isalpha(temp[0]) ) {
			continue;
		}
		middle = temp.find(" ");
		padding.insert(std::make_pair(temp.substr(0, middle), isa::utils::castToType< std::string, unsigned int >(temp.substr(middle + 1))));
	}
}

void readVectorWidth(map< std::string, unsigned int > & vectorWidth, const std::string & vectorFilename) {
	std::string temp;
	std::ifstream vectorFile(vectorFilename);

	while ( ! vectorFile.eof() ) {
		unsigned int middle = 0;

		std::getline(vectorFile, temp);
		if ( ! std::isalpha(temp[0]) ) {
			continue;
		}
		middle = temp.find(" ");
		vectorWidth.insert(std::make_pair(temp.substr(0, middle), isa::utils::castToType< std::string, unsigned int >(temp.substr(middle + 1))));
	}
}

void readDedispersion(map< std::string, map< unsigned int, std::vector< unsigned int > > > & dedispersionParameters, const std::string  & dedispersionFilename) {
	std::string temp;
	std::ifstream dedispersionFile(dedispersionFilename);

	while ( ! dedispersionFile.eof() ) {
		unsigned int splitPoint = 0;

		std::getline(dedispersionFile, temp);
		if ( ! std::isalpha(temp[0]) ) {
			continue;
		}
		std::string deviceName;
		unsigned int nrDMs = 0;
		std::vector< unsigned int > parameters(4);

		splitPoint = temp.find(" ");
		deviceName = temp.substr(0, splitPoint);
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		nrDMs = isa::utils::castToType< std::string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		parameters[0] = isa::utils::castToType< std::string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		parameters[1] = isa::utils::castToType< std::string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		parameters[2] = isa::utils::castToType< std::string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		parameters[3] = isa::utils::castToType< std::string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		parameters[4] = isa::utils::castToType< std::string, unsigned int >(temp);

		if ( dedispersionParameters.count(deviceName) == 0 ) {
			map< unsigned int, std::vector< unsigned int > > container;

			container.insert(std::make_pair(nrDMs, parameters));
			dedispersionParameters.insert(std::make_pair(deviceName, container));
		} else {
			dedispersionParameters[deviceName].insert(std::make_pair(nrDMs, parameters));
		}
	}
}

void readTranspose(map< std::string, map< unsigned int, unsigned int > > & transposeParameters, const std::string & transposeFilename) {
	std::string temp;
	std::ifstream transposeFile(transposeFilename);

	while ( ! transposeFile.eof() ) {
		unsigned int splitPoint = 0;

		std::getline(transposeFile, temp);
		if ( ! std::isalpha(temp[0]) ) {
			continue;
		}
		std::string deviceName;
		unsigned int nrDMs = 0;
		unsigned int parameter = 0;

		splitPoint = temp.find(" ");
		deviceName = temp.substr(0, splitPoint);
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		nrDMs = isa::utils::castToType< std::string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		parameter = isa::utils::castToType< std::string, unsigned int >(temp);

		if ( transposeParameters.count(deviceName) == 0 ) {
			map< unsigned int, unsigned int > container;

			container.insert(std::make_pair(nrDMs, parameter));
			transposeParameters.insert(std::make_pair(deviceName, container));
		} else {
			transposeParameters[deviceName].insert(std::make_pair(nrDMs, parameter));
		}
	}
}

void readFolding(map< std::string, map< unsigned int, map< unsigned int, std::vector< unsigned int > > > > & foldingParameters, const std::string & foldingFilename) {
	std::string temp;
	std::ifstream foldingFile(foldingFilename);

	while ( ! foldingFile.eof() ) {
		unsigned int splitPoint = 0;

		std::getline(foldingFile, temp);
		if ( ! std::isalpha(temp[0]) ) {
			continue;
		}
		std::string deviceName;
		unsigned int nrDMs = 0;
		unsigned int nrPeriods = 0;
		std::vector< unsigned int > parameters(6);

		splitPoint = temp.find(" ");
		deviceName = temp.substr(0, splitPoint);
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		nrDMs = isa::utils::castToType< std::string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		nrPeriods = isa::utils::castToType< std::string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		parameters[0] = isa::utils::castToType< std::string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		parameters[1] = isa::utils::castToType< std::string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		parameters[2] = isa::utils::castToType< std::string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		parameters[3] = isa::utils::castToType< std::string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		parameters[4] = isa::utils::castToType< std::string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		parameters[5] = isa::utils::castToType< std::string, unsigned int >(temp);

		if ( foldingParameters.count(deviceName) == 0 ) {
			map< unsigned int, map< unsigned int, std::vector< unsigned int > > > externalContainer;
			map< unsigned int, std::vector< unsigned int > > internalContainer;

			internalContainer.insert(std::make_pair(nrPeriods, parameters));
			externalContainer.insert(std::make_pair(nrDMs, internalContainer));
			foldingParameters.insert(std::make_pair(deviceName, externalContainer));
		} else if ( foldingParameters[deviceName].count(nrDMs) == 0 ) {
			map< unsigned int, std::vector< unsigned int > > internalContainer;

			internalContainer.insert(std::make_pair(nrPeriods, parameters));
			foldingParameters[deviceName].insert(std::make_pair(nrDMs, internalContainer));
		} else {
			foldingParameters[deviceName][nrDMs].insert(std::make_pair(nrPeriods, parameters));
		}
	}
}

void readSNR(map< std::string, map< unsigned int, map< unsigned int, std::vector< unsigned int > > > > & snrParameters, const std::string & snrFilename) {
	std::string temp;
	std::ifstream snrFile(snrFilename);

	while ( ! snrFile.eof() ) {
		unsigned int splitPoint = 0;

		std::getline(snrFile, temp);
		if ( ! std::isalpha(temp[0]) ) {
			continue;
		}
		std::string deviceName;
		unsigned int nrDMs = 0;
		unsigned int nrPeriods = 0;
		std::vector< unsigned int > parameters(2);

		splitPoint = temp.find(" ");
		deviceName = temp.substr(0, splitPoint);
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		nrDMs = isa::utils::castToType< std::string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		nrPeriods = isa::utils::castToType< std::string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		parameters[0] = isa::utils::castToType< std::string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		parameters[1] = isa::utils::castToType< std::string, unsigned int >(temp);

		if ( snrParameters.count(deviceName) == 0 ) {
			map< unsigned int, map< unsigned int, std::vector< unsigned int > > > externalContainer;
			map< unsigned int, std::vector< unsigned int > > internalContainer;

			internalContainer.insert(std::make_pair(nrPeriods, parameters));
			externalContainer.insert(std::make_pair(nrDMs, internalContainer));
			snrParameters.insert(std::make_pair(deviceName, externalContainer));
		} else if ( snrParameters[deviceName].count(nrDMs) == 0 ) {
			map< unsigned int, std::vector< unsigned int > > internalContainer;

			internalContainer.insert(std::make_pair(nrPeriods, parameters));
			snrParameters[deviceName].insert(std::make_pair(nrDMs, internalContainer));
		} else {
			snrParameters[deviceName][nrDMs].insert(std::make_pair(nrPeriods, parameters));
		}
	}
}

