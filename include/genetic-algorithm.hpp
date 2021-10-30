#pragma once
///@file genetic-algorithm.hpp
///@brief Library that provides a variety of genetic algorithms
#include <cstdint>

/*!
 * @brief Picks winners based on a roulette wheel whose sectors' widths are proportional to the fitness of each individual
 * @param[in]   populationSize  The number of individuals
 * @param[in]   fitness         Array to the fitnesses of each individual
 * @param[out]  winners         Array that will be filled with the indices of the picked winners
 * @param[in]   winnersSize     The desired number of winners
 */
void rouletteRanking(int populationSize, float *fitness, int *winners, int winnersSize);

/*!
 * @brief Picks winners based on a roulette wheel whose sectors' width is based on the ranks of the individuals, with selectionPressure acting as a parameter to determine how much rank weighs.
 * The formula used is P(r) = k1 - r*k2 where r=rank, k1= selectionPressure/populationSize and k2=selectionPressure/(populationSize*(populationSize-1))
 * @param[in]   populationSize      The number of individuals
 * @param[in]   fitness             Array to the fitnesses of each individual
 * @param[in]   maximizeFitness     True if the objective is to maximize fitness, False otherwise
 * @param[in]   selectionPressure   Determines how much each rank weighs
 * @param[out]  winners             Array that will be filled with the indices of the picked winners
 * @param[in]   winnersSize         The desired number of winners
 */
void linearRanking(int populationSize, float *fitness, bool maximizeFitness, float selectionPressure, int *winners, int winnersSize);

/*!
 * @brief Picks winners based on a roulette wheel whose sectors' width is based on the ranks of the individuals, with k1 acting as a parameter to determine how much rank weighs.
 * The formula used is P(r) = k1*(1-k1)^r where r=rank
 * @param[in]   populationSize      The number of individuals
 * @param[in]   fitness             Array to the fitnesses of each individual
 * @param[in]   maximizeFitness     True if the objective is to maximize fitness, False otherwise
 * @param[in]   k1                  Determines how much each rank weighs
 * @param[out]  winners             Array that will be filled with the indices of the picked winners
 * @param[in]   winnersSize         The desired number of winners
 */
void exponentialRanking(int populationSize, float *fitness, bool maximizeFitness, float k1, int *winners, int winnersSize);

/*!
 * @brief Picks winners by randomly selecting tournamentSize individuals winnersSize times, and then picking the individual with the highest fitness in each tournament
 * @param[in]   populationSize      The number of individuals
 * @param[in]   fitness             Array to the fitnesses of each individual
 * @param[in]   maximizeFitness     True if the objective is to maximize fitness, False otherwise
 * @param[in]   tournamentSize      Determines the size of each tournament
 * @param[out]  winners             Array that will be filled with the indices of the picked winners
 * @param[in]   winnersSize         The desired number of winners
 */
void tournamentRanking(int populationSize, float *fitness, bool maximizeFitness, int tournamentSize, int *winners, int winnersSize);

/*!
 * @brief Selects two random points among the ones defined in genesLoci and uses them to cut up and paste together three alternating sections from the two parents. Use this if not all of your genes are 1 byte long
 * @param[in]   parent1         The first parent
 * @param[in]   parent2         The second parent
 * @param[in]   length          The length of the parents and, by consequence, the child
 * @param[out]  child           The result of cutting up and pasting the parents
 * @param[in]   genesLoci       The starting and ending positions of each gene
 * @note While the function has code to recognize whether you've included the extremes (0, length-1), it would need to allocate memory for them, so it's better if you generate genesLoci as having them yourself
 * @param[in]   genesLociLength The length of genesLoci
 */
void twoPointsCrossover(uint8_t *parent1, uint8_t *parent2, uint64_t length, uint8_t *child, uint64_t *genesLoci, int genesLociLength);

/*!
 * @brief Selects two random points and uses them to cut up and paste together three alternating sections from the two parents. Use this if every gene is 1 byte long.
 * @param[in]   parent1         The first parent
 * @param[in]   parent2         The second parent
 * @param[in]   length          The length of the parents and, by consequence, the child
 * @param[out]  child           The result of cutting up and pasting the parents
 */
void twoPointsCrossover(uint8_t *parent1, uint8_t *parent2, uint64_t length, uint8_t *child);

/*!
 * @brief For each gene, as defined by genesLoci, selects whether child will inherit it from parent1 or parent2. Use this if not all of your genes are 1 byte long
 * @param[in]   parent1         The first parent
 * @param[in]   parent2         The second parent
 * @param[in]   length          The length of the parents and, by consequence, the child
 * @param[out]  child           The result of cutting up and pasting the parents
 * @param[in]   genesLoci       The starting and ending positions of each gene
 * @note While the function has code to recognize whether you've included the extremes (0, length-1), it would need to allocate memory for them, so it's better if you generate genesLoci as having them yourself
 * @param[in]   genesLociLength The length of genesLoci
 */
void uniformCrossover(uint8_t *parent1, uint8_t *parent2, uint64_t length, uint8_t *child, uint64_t *genesLoci, int genesLociLength);


/*!
 * @brief For each gene selects whether child will inherit it from parent1 or parent2. Use this if all of your genes are 1 byte long
 * @param[in]   parent1         The first parent
 * @param[in]   parent2         The second parent
 * @param[in]   length          The length of the parents and, by consequence, the child
 * @param[out]  child           The result of cutting up and pasting the parents
 */
void uniformCrossover(uint8_t *parent1, uint8_t *parent2, uint64_t length, uint8_t *child);

/*!
 * @brief Alters a random bit of each byte with mutationProbability probability. Use this if your genome has only 1 byte long genes, and all possible values for the genes are accepted. Otherwise you'll have to define your own function
 * @param[out]  individual          The genome to mutate
 * @param[in]   length              The length of individual
 * @param[in]   mutationProbability The probability each gene mutates
 */
void mutate(uint8_t *individual, int length, float mutationProbability);
