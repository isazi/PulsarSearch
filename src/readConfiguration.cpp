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
using std::getline;
#include <fstream>
using std::ifstream;
#include <map>
using std::make_pair;
#include <vector>
using std::vector;
#include <cctype>
using std::isalpha;

#include <readConfiguration.hpp>
#include <utils.hpp>
using isa::utils::castToType;


void readPadding(map< string, unsigned int > & padding, const string paddingFilename) {
	string temp;
	ifstream paddingFile(paddingFilename);

	while ( ! paddingFile.eof() ) {
		unsigned int middle = 0;

		getline(paddingFile, temp);
		if ( ! isalpha(temp[0]) ) {
			continue;
		}
		middle = temp.find(" ");
		padding.insert(make_pair(temp.substr(0, middle), castToType< string, unsigned int >(temp.substr(middle + 1))));
	}
}

void readVectorWidth(map< string, unsigned int > & vectorWidth, const string vectorFilename) {
	string temp;
	ifstream vectorFile(vectorFilename);

	while ( ! vectorFile.eof() ) {
		unsigned int middle = 0;

		getline(vectorFile, temp);
		if ( ! isalpha(temp[0]) ) {
			continue;
		}
		middle = temp.find(" ");
		vectorWidth.insert(make_pair(temp.substr(0, middle), castToType< string, unsigned int >(temp.substr(middle + 1))));
	}
}

void readDedispersion(map< string, map< unsigned int, vector< unsigned int > > > & dedispersionParameters, const string dedispersionFilename) {
	string temp;
	ifstream dedispersionFile(dedispersionFilename);

	while ( ! dedispersionFile.eof() ) {
		unsigned int splitPoint = 0;

		getline(dedispersionFile, temp);
		if ( ! isalpha(temp[0]) ) {
			continue;
		}
		string deviceName;
		unsigned int nrDMs = 0;
		vector< unsigned int > parameters(4);

		splitPoint = temp.find(" ");
		deviceName = temp.substr(0, splitPoint);
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		nrDMs = castToType< string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		parameters[0] = castToType< string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		parameters[1] = castToType< string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		parameters[2] = castToType< string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		parameters[3] = castToType< string, unsigned int >(temp);

		if ( dedispersionParameters.count(deviceName) == 0 ) {
			map< unsigned int, vector< unsigned int > > container;

			container.insert(make_pair(nrDMs, parameters));
			dedispersionParameters.insert(make_pair(deviceName, container));
		} else {
			dedispersionParameters[deviceName].insert(make_pair(nrDMs, parameters));
		}
	}
}

void readTranspose(map< string, map< unsigned int, unsigned int > > & transposeParameters, const string transposeFilename) {
	string temp;
	ifstream transposeFile(transposeFilename);

	while ( ! transposeFile.eof() ) {
		unsigned int splitPoint = 0;

		getline(transposeFile, temp);
		if ( ! isalpha(temp[0]) ) {
			continue;
		}
		string deviceName;
		unsigned int nrDMs = 0;
		unsigned int parameter = 0;

		splitPoint = temp.find(" ");
		deviceName = temp.substr(0, splitPoint);
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		nrDMs = castToType< string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		parameter = castToType< string, unsigned int >(temp);

		if ( transposeParameters.count(deviceName) == 0 ) {
			map< unsigned int, unsigned int > container;

			container.insert(make_pair(nrDMs, parameter));
			transposeParameters.insert(make_pair(deviceName, container));
		} else {
			transposeParameters[deviceName].insert(make_pair(nrDMs, parameter));
		}
	}
}

void readFolding(map< string, map< unsigned int, map< unsigned int, vector< unsigned int > > > > & foldingParameters, const string foldingFilename) {
	string temp;
	ifstream foldingFile(foldingFilename);

	while ( ! foldingFile.eof() ) {
		unsigned int splitPoint = 0;

		getline(foldingFile, temp);
		if ( ! isalpha(temp[0]) ) {
			continue;
		}
		string deviceName;
		unsigned int nrDMs = 0;
		unsigned int nrPeriods = 0;
		vector< unsigned int > parameters(6);

		splitPoint = temp.find(" ");
		deviceName = temp.substr(0, splitPoint);
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		nrDMs = castToType< string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		nrPeriods = castToType< string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		parameters[0] = castToType< string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		parameters[1] = castToType< string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		parameters[2] = castToType< string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		parameters[3] = castToType< string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		parameters[4] = castToType< string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		parameters[5] = castToType< string, unsigned int >(temp);

		if ( foldingParameters.count(deviceName) == 0 ) {
			map< unsigned int, map< unsigned int, vector< unsigned int > > > externalContainer;
			map< unsigned int, vector< unsigned int > > internalContainer;

			internalContainer.insert(make_pair(nrPeriods, parameters));
			externalContainer.insert(make_pair(nrDMs, internalContainer));
			foldingParameters.insert(make_pair(deviceName, externalContainer));
		} else if ( foldingParameters[deviceName].count(nrDMs) == 0 ) {
			map< unsigned int, vector< unsigned int > > internalContainer;

			internalContainer.insert(make_pair(nrPeriods, parameters));
			foldingParameters[deviceName].insert(make_pair(nrDMs, internalContainer));
		} else {
			foldingParameters[deviceName][nrDMs].insert(make_pair(nrPeriods, parameters));
		}
	}
}

void readSNR(map< string, map< unsigned int, map< unsigned int, vector< unsigned int > > > > & snrParameters, const string snrFilename) {
	string temp;
	ifstream snrFile(snrFilename);

	while ( ! snrFile.eof() ) {
		unsigned int splitPoint = 0;

		getline(snrFile, temp);
		if ( ! isalpha(temp[0]) ) {
			continue;
		}
		string deviceName;
		unsigned int nrDMs = 0;
		unsigned int nrPeriods = 0;
		vector< unsigned int > parameters(2);

		splitPoint = temp.find(" ");
		deviceName = temp.substr(0, splitPoint);
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		nrDMs = castToType< string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		nrPeriods = castToType< string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		parameters[0] = castToType< string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		parameters[1] = castToType< string, unsigned int >(temp);

		if ( snrParameters.count(deviceName) == 0 ) {
			map< unsigned int, map< unsigned int, vector< unsigned int > > > externalContainer;
			map< unsigned int, vector< unsigned int > > internalContainer;

			internalContainer.insert(make_pair(nrPeriods, parameters));
			externalContainer.insert(make_pair(nrDMs, internalContainer));
			snrParameters.insert(make_pair(deviceName, externalContainer));
		} else if ( snrParameters[deviceName].count(nrDMs) == 0 ) {
			map< unsigned int, vector< unsigned int > > internalContainer;

			internalContainer.insert(make_pair(nrPeriods, parameters));
			snrParameters[deviceName].insert(make_pair(nrDMs, internalContainer));
		} else {
			snrParameters[deviceName][nrDMs].insert(make_pair(nrPeriods, parameters));
		}
	}
}
