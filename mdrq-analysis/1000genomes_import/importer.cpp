/*
 * Importer for transforming a 1000 Genomes Project VCF file into features vectors.
 * VCF files can be downloaded from: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
 *
 * Author: Stefan Sprenger
 */

#define DIM 19
#define SAMPLES 2504
#define PED_FILE "samples.ped"
#define REFERENCE_GENOME "hs37d5"
#define VECTORS 10000000

#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

struct Sample {
  std::string id, population, family_id, relationship;
  uint32_t gender;
} Sample;

struct Mutation {
  uint32_t chromosome, position, depth, quality, allele_count;
  std::string variant_id, genotype, reference_genome, filter,
              ancestral_allele, variant_type, ref_base, alt_base;
  float allele_freq;
} Mutation;

/**
 * Returns a order-preserving hash value.
 * Extracts a digit from the source value and combines it with
 * first n values of its hash value (n is set to 5 by default).
 *
 * For instance, the string 'HG1013' becomes: 101332145
 */
std::size_t getOrderPreservingHashValue(std::string sourceValue) {
  std::stringstream stringStream;
  std::size_t orderPreservingHashValue;
  std::string formatColumn = sourceValue;
  std::smatch match;

  // get digit from source value (max 10000)
  std::regex rgx("[0-9]+");
  if (std::regex_search(formatColumn, match, rgx))
    stringStream << std::to_string((stoi(match[0].str()) % 10000));

  // get first 5 values of hash of source value
  std::size_t str_hash  = std::hash<std::string>()(sourceValue);
  std::string hashValue = std::to_string(str_hash);
  stringStream << hashValue.substr(0,4);

  // transform back into integer
  orderPreservingHashValue = stoi(stringStream.str());

  return orderPreservingHashValue;
}

/**
 * Transforms a combination of a mutation and a sample record
 * into a feature vector which is basically an array of DIM floats.
 */
std::vector<float> getFeatureVector(struct Mutation mutation, struct Sample sample) {
  std::vector<float> featureVector(DIM);
  std::size_t str_hash;

  featureVector[0]  = (float) getOrderPreservingHashValue(sample.id);
  str_hash = std::hash<std::string>()(sample.population);
  featureVector[1]  = (float) str_hash;
  featureVector[2]  = (float) sample.gender;
  featureVector[3]  = (float) getOrderPreservingHashValue(sample.family_id);
  str_hash = std::hash<std::string>()(sample.relationship);
  featureVector[4]  = (float) str_hash;
  featureVector[5]  = (float) mutation.chromosome;
  featureVector[6]  = (float) mutation.position;
  featureVector[7]  = (float) getOrderPreservingHashValue(mutation.variant_id);
  featureVector[8]  = (float) mutation.quality;
  featureVector[9]  = (float) mutation.depth;
  featureVector[10] = mutation.allele_freq;
  featureVector[11] = (float) mutation.allele_count;
  str_hash = std::hash<std::string>()(mutation.filter);
  featureVector[12] = (float) str_hash;
  str_hash = std::hash<std::string>()(mutation.ref_base);
  featureVector[13] = (float) str_hash;
  str_hash = std::hash<std::string>()(mutation.alt_base);
  featureVector[14] = (float) str_hash;
  str_hash = std::hash<std::string>()(mutation.variant_type);
  featureVector[15] = (float) str_hash;
  str_hash = std::hash<std::string>()(mutation.ancestral_allele);
  featureVector[16] = (float) str_hash;
  str_hash = std::hash<std::string>()(mutation.genotype);
  featureVector[17] = (float) str_hash;
  str_hash = std::hash<std::string>()(mutation.reference_genome);
  featureVector[18] = (float) str_hash;

  return featureVector;
}

void process_sample_names(struct Sample *samples, std::vector<std::string> line_tokens) {
  for (size_t i = 9; i < line_tokens.size(); ++i) {
    size_t j = i-9;
    samples[j].id = line_tokens[i];
    std::ifstream peds(PED_FILE);
    std::string line;
    std::string token;

    while (std::getline(peds, line)) {
      std::vector<std::string> line_tokens;
      std::istringstream iss(line);
      while(std::getline(iss, token, '\t'))
        line_tokens.push_back(token);
      if (line_tokens[1] == samples[j].id) {
        samples[j].gender       = stoi(line_tokens[4]);
        samples[j].population   = line_tokens[6];
        samples[j].family_id    = line_tokens[0];
        samples[j].relationship = line_tokens[7];
      }
    }
  }
  std::cout << "Imported sample metadata from " << PED_FILE << "." << std::endl;
}

int main(int argc, char* argv[]) {
  if (argc != 3) {
    std::cout << "Usage: " << argv[0] << " vcf_filename featurevector_filename"
              << std::endl;
    return 1;
  }

  // variables needed for file processing
  size_t i = 0;
  std::string line;
  std::string vcf_file = argv[1];
  std::string fvector_file = argv[2];
  std::ifstream input_stream(vcf_file);
  std::ofstream output_stream(fvector_file);

  struct Sample *samples = new struct Sample[SAMPLES];
  std::cout << "Starting transformation of '" << vcf_file << "' into '"
            << fvector_file << "'." << std::endl;

  if (input_stream) {
    std::string token;
    // Linewise processing of VCF file
    while (std::getline(input_stream, line) && i++ < VECTORS) {
      std::vector<std::string> line_tokens;
      struct Mutation mutation;
      std::istringstream iss(line);
      // VCF uses tabs as delimiter
      while(std::getline(iss, token, '\t'))
        line_tokens.push_back(token);

      // ignore comments
      if (line_tokens[0][0] == '#' && line_tokens[0][1] == '#')
        continue;
      // process sample names
      if (line_tokens[0][0] == '#') {
        process_sample_names(samples, line_tokens);
        continue;
      }

      mutation.chromosome = stoi(line_tokens[0]);
      mutation.position   = stoi(line_tokens[1]);
      mutation.variant_id = line_tokens[2];
      mutation.ref_base   = line_tokens[3];
      mutation.alt_base   = line_tokens[4];
      mutation.quality    = stoi(line_tokens[5]);
      mutation.filter     = line_tokens[6];
      mutation.reference_genome = REFERENCE_GENOME;

      // process FORMAT column using regular expressions
      std::string formatColumn = line_tokens[7];
      std::smatch match;

      // parse allele frequency
      std::regex rgx("AF=0.[0-9]+");
      if (std::regex_search(formatColumn, match, rgx))
        mutation.allele_freq = stof(match[0].str().substr(3));

      // parse depth
      rgx = std::regex("DP=[0-9]+");
      if (std::regex_search(formatColumn, match, rgx))
        mutation.depth = stoi(match[0].str().substr(3));

      // parse ancestral allele
      rgx = std::regex("AA=[.|]+");
      if (std::regex_search(formatColumn, match, rgx))
        mutation.ancestral_allele = match[0].str().substr(3);

      // parse variant_type
      rgx = std::regex("VT=[a-zA-Z]+");
      if (std::regex_search(formatColumn, match, rgx))
        mutation.variant_type = match[0].str().substr(3);

      // parse allele_count
      rgx = std::regex("AC=[0-9]+");
      if (std::regex_search(formatColumn, match, rgx))
        mutation.allele_count = stoi(match[0].str().substr(3));

      // process samples
      for (size_t j = 9; j < line_tokens.size(); ++j) {
        // ignore samples without the mutation
        if (line_tokens[j] == "0|0")
          continue;

        mutation.genotype = line_tokens[j];

        std::vector<float> featureVector = getFeatureVector(mutation, samples[j-9]);
        // write feature vector into output file
        output_stream << featureVector[0];
        for (size_t i = 1; i < featureVector.size(); ++i)
          output_stream << " " << featureVector[i];
        output_stream << std::endl;
      }
      if (i % 1000 == 0)
	std::cout << "Imported " << i << " vectors." << std::endl;
    }

    std::cout << std::endl << "SUCCESS: File '" << vcf_file
              << "' has been transformed into features vectors." << std::endl;
  } else {
    std::cout << std::endl << "ERROR: File '" << vcf_file
              << "' could not be found." << std::endl;
  }

  input_stream.close();
  output_stream.close();

  return 0;
}
