
#include <Dedispersion.hpp>

#ifndef CONFIGURATION_HPP
#define CONFIGURATION_HPP

// Types for the data
typedef float dataType;
const std::string dataName("float");

// DEBUG mode, prints to screen some useful information
const bool DEBUG = true;

// Memory padding
std::map< std::string, unsigned int > padding;

// Vector unit width
std::map< std::string, unsigned int > vectorWidth;

// Tuned parameters
std::map< std::string, std::map< unsigned int, PulsarSearch::DedispersionConf > > dedispersionParameters;
std::map< std::string, std::map< unsigned int, isa::OpenCL::transposeConf > > transposeParameters;
std::map< std::string, std::map< unsigned int, PulsarSearch::snrDedispersedConf > > snrDParameters;
std::map< std::string, std::map< unsigned int, std::map< unsigned int, PulsarSearch::FoldingConf > > > foldingParameters;
std::map< std::string, std::map< unsigned int, std::map< unsigned int, PulsarSearch::snrFoldedConf > > > snrFDParameters;

#endif // CONFIGURATION_HPP

