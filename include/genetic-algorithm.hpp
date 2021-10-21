#pragma once
#include <cstdint>

void rouletteRanking(int populationSize, float *fitness, int *winners, int winnersSize);

void linearRanking(int populationSize, float *fitness, bool maximizeFitness, float selectionPressure, int *winners, int winnersSize);

void exponentialRanking(int populationSize, float *fitness, bool maximizeFitness, float k1, int *winners, int winnersSize);

void tournamentRanking(int populationSize, float *fitness, bool maximizeFitness, int tournamentSize, int *winners, int winnersSize);

void twoPointsCrossover(uint8_t *parent1, uint8_t *parent2, uint64_t length, uint8_t *child, uint64_t *genesLoci, int genesLociLength);

void twoPointsCrossover(uint8_t *parent1, uint8_t *parent2, uint64_t length, uint8_t *child);

void uniformCrossover(uint8_t *parent1, uint8_t *parent2, uint64_t length, uint8_t *child);

void uniformCrossover(uint8_t *parent1, uint8_t *parent2, uint64_t length, uint8_t *child, uint64_t *genesLoci, int genesLociLength);

void mutate(uint8_t *individual, int length, float mutationProbability);
