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


void readPadding(std::map< std::string, unsigned int > & padding, const std::string & paddingFilename) {
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

void readVectorWidth(std::map< std::string, unsigned int > & vectorWidth, const std::string & vectorFilename) {
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

void readDedispersion(std::map< std::string, std::map< unsigned int, PulsarSearch::DedispersionConf > > & dedispersionParameters, const std::string  & dedispersionFilename) {
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
    PulsarSearch::DedispersionConf conf;

		splitPoint = temp.find(" ");
		deviceName = temp.substr(0, splitPoint);
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		nrDMs = isa::utils::castToType< std::string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		conf.setLocalMem(isa::utils::castToType< std::string, bool >(temp.substr(0, splitPoint)));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		conf.setUnroll(isa::utils::castToType< std::string, unsigned int >(temp.substr(0, splitPoint)));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		conf.setNrSamplesPerBlock(isa::utils::castToType< std::string, unsigned int >(temp.substr(0, splitPoint)));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		conf.setNrDMsPerBlock(isa::utils::castToType< std::string, unsigned int >(temp.substr(0, splitPoint)));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		conf.setNrSamplesPerThread(isa::utils::castToType< std::string, unsigned int >(temp.substr(0, splitPoint)));
		temp = temp.substr(splitPoint + 1);
		conf.setNrDMsPerThread(isa::utils::castToType< std::string, unsigned int >(temp));

		if ( dedispersionParameters.count(deviceName) == 0 ) {
      std::map< unsigned int, PulsarSearch::DedispersionConf > container;

			container.insert(std::make_pair(nrDMs, conf));
			dedispersionParameters.insert(std::make_pair(deviceName, container));
		} else {
			dedispersionParameters[deviceName].insert(std::make_pair(nrDMs, conf));
		}
	}
}

void readTranspose(std::map< std::string, std::map< unsigned int, isa::OpenCL::transposeConf > > & transposeParameters, const std::string & transposeFilename) {
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
    isa::OpenCL::transposeConf parameters;

		splitPoint = temp.find(" ");
		deviceName = temp.substr(0, splitPoint);
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		nrDMs = isa::utils::castToType< std::string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		parameters.setNrItemsPerBlock(isa::utils::castToType< std::string, unsigned int >(temp));

		if ( transposeParameters.count(deviceName) == 0 ) {
      std::map< unsigned int, isa::OpenCL::transposeConf > container;

			container.insert(std::make_pair(nrDMs, parameters));
			transposeParameters.insert(std::make_pair(deviceName, container));
		} else {
			transposeParameters[deviceName].insert(std::make_pair(nrDMs, parameters));
		}
	}
}

void readSNRD(std::map< std::string, std::map< unsigned int, PulsarSearch::snrDedispersedConf > > & snrParameters, const std::string & snrFilename) {
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
    PulsarSearch::snrDedispersedConf parameters;

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
		parameters.setNrDMsPerBlock(isa::utils::castToType< std::string, unsigned int >(temp.substr(0, splitPoint)));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		parameters.setNrDMsPerThread(isa::utils::castToType< std::string, unsigned int >(temp.substr(0, splitPoint)));

		if ( snrParameters.count(deviceName) == 0 ) {
      std::map< unsigned int, PulsarSearch::snrDedispersedConf > container;

			container.insert(std::make_pair(nrDMs, parameters));
			snrParameters.insert(std::make_pair(deviceName, container));
		} else {
			snrParameters[deviceName].insert(std::make_pair(nrDMs, parameters));
		}
	}
}

void readFolding(std::map< std::string, std::map< unsigned int, std::map< unsigned int, PulsarSearch::FoldingConf > > > & foldingParameters, const std::string & foldingFilename) {
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
    PulsarSearch::FoldingConf parameters;

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
		parameters.setNrDMsPerBlock(isa::utils::castToType< std::string, unsigned int >(temp.substr(0, splitPoint)));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		parameters.setNrPeriodsPerBlock(isa::utils::castToType< std::string, unsigned int >(temp.substr(0, splitPoint)));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		parameters.setNrBinsPerBlock(isa::utils::castToType< std::string, unsigned int >(temp.substr(0, splitPoint)));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		parameters.setNrDMsPerThread(isa::utils::castToType< std::string, unsigned int >(temp.substr(0, splitPoint)));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		parameters.setNrPeriodsPerThread(isa::utils::castToType< std::string, unsigned int >(temp.substr(0, splitPoint)));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		parameters.setNrBinsPerThread(isa::utils::castToType< std::string, unsigned int >(temp.substr(0, splitPoint)));
		temp = temp.substr(splitPoint + 1);
		parameters.setVector(isa::utils::castToType< std::string, unsigned int >(temp));

		if ( foldingParameters.count(deviceName) == 0 ) {
      std::map< unsigned int, std::map< unsigned int, PulsarSearch::FoldingConf > > externalContainer;
      std::map< unsigned int, PulsarSearch::FoldingConf > internalContainer;

			internalContainer.insert(std::make_pair(nrPeriods, parameters));
			externalContainer.insert(std::make_pair(nrDMs, internalContainer));
			foldingParameters.insert(std::make_pair(deviceName, externalContainer));
		} else if ( foldingParameters[deviceName].count(nrDMs) == 0 ) {
      std::map< unsigned int, std::vector< unsigned int > > internalContainer;

			internalContainer.insert(std::make_pair(nrPeriods, parameters));
			foldingParameters[deviceName].insert(std::make_pair(nrDMs, internalContainer));
		} else {
			foldingParameters[deviceName][nrDMs].insert(std::make_pair(nrPeriods, parameters));
		}
	}
}

void readSNRF(std::map< std::string, std::map< unsigned int, std::map< unsigned int, PulsarSearch::snrFoldedConf > > > & snrParameters, const std::string & snrFilename) {
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
    PulsarSearch::snrFoldedConf parameters;

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
		parameters.setNrDMsPerBlock(isa::utils::castToType< std::string, unsigned int >(temp.substr(0, splitPoint)));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		parameters.setNrPeriodsPerBlock(isa::utils::castToType< std::string, unsigned int >(temp.substr(0, splitPoint)));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		parameters.setNrDMsPerThread(isa::utils::castToType< std::string, unsigned int >(temp.substr(0, splitPoint)));
		temp = temp.substr(splitPoint + 1);
		parameters.setNrPeriodsPerThread(isa::utils::castToType< std::string, unsigned int >(temp));

		if ( snrParameters.count(deviceName) == 0 ) {
      std::map< unsigned int, std::map< unsigned int, PulsarSearch::snrFoldedConf > > externalContainer;
      std::map< unsigned int, PulsarSearch::snrFoldedConf > internalContainer;

			internalContainer.insert(std::make_pair(nrPeriods, parameters));
			externalContainer.insert(std::make_pair(nrDMs, internalContainer));
			snrParameters.insert(std::make_pair(deviceName, externalContainer));
		} else if ( snrParameters[deviceName].count(nrDMs) == 0 ) {
      std::map< unsigned int, std::vector< unsigned int > > internalContainer;

			internalContainer.insert(std::make_pair(nrPeriods, parameters));
			snrParameters[deviceName].insert(std::make_pair(nrDMs, internalContainer));
		} else {
			snrParameters[deviceName][nrDMs].insert(std::make_pair(nrPeriods, parameters));
		}
	}
}

