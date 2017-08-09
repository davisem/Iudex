/*
 * fastq_reader.cpp
 *
 *  Created on: Aug 8, 2017
 *      Author: ericdavis
 */

#include <iostream>
#include <fstream>
#include <istream>
#include "bloom_filter.hpp"

using namespace std;

int main (int argc, char* argv[] ) {
	int filtered_count = 0;

	string name;
	string seq;
	string plus;
	string quals;
	string line;

	ifstream input_fastq(argv[1]);
	ofstream filtered_fastq(argv[2]);
	double false_positive_probability = atof (argv[3]);

	bloom_parameters parameters;
	parameters.projected_element_count = 10000000;
	parameters.false_positive_probability = false_positive_probability;
	parameters.random_seed = 0xA5A5A5A5;

	parameters.compute_optimal_parameters();

	bloom_filter filter(parameters);

	while (std::getline(input_fastq, line)) {

		name = line;

		if (name.at(0) !='@') {
			cout << "Malformed fastq file, expected @, found: " << line << endl;
			return 1;
			}

		getline(input_fastq, seq);
		getline(input_fastq, plus);
		getline(input_fastq, quals);

		if (filter.contains(seq)) {
			++ filtered_count;
			continue;
		}

		filter.insert(seq);

		filtered_fastq << name << endl;
		filtered_fastq << seq << endl;
		filtered_fastq << plus << endl;
		filtered_fastq << quals << endl;

	}

	cout << "Number of filtered records:" << filtered_count << endl;
	return 0;
}
