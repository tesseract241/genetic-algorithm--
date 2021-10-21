#include <cassert>
#include <cstdint>
#include <algorithm>
#include <random>
#include <cstring>
#include <unordered_map>
#include <vector>
#include <genetic-algorithm.hpp>

const int storedProbabilities   = 5;
const int trackedProbabilities  = 10;
struct RanksData{
    int populationSize;
    int usage;
    float* probabilities;
};

std::unordered_map<float, std::unordered_map<int, RanksData>> linearRanksData;
std::unordered_map<float, RanksData> exponentialRanksData;

void rouletteRanking(int populationSize, float *fitness, int *winners, int winnersSize){
    assert(populationSize && "rouletteRanking: populationSize was 0\n");
    float *cumulativeProbabilities = new float[populationSize];
    float minFitness = *std::min_element(fitness, fitness+populationSize);
    if(minFitness>=0){
        cumulativeProbabilities[0] = fitness[0];
        for(int i=1;i<populationSize;++i){
            cumulativeProbabilities[i] = cumulativeProbabilities[i-1] + fitness[i];
        }
        float reciprocal_sum = 1./fitness[populationSize-1];
        for (int i=0;i<populationSize;++i){
            cumulativeProbabilities[i] *= reciprocal_sum;
        }
    } else{
        //For problems in which the fitness must be minimized, the modified fitness = 1/(1 + fitness - minFitness(len(population))
        //is used.
        minFitness/=float(populationSize);
        cumulativeProbabilities[0] = 1./(1. + fitness[0] - minFitness);
        for(int i=1;i<populationSize;++i){
            cumulativeProbabilities[i] = cumulativeProbabilities[i-1] + 1./(1. + fitness[i] - minFitness);
        }
        float reciprocal_sum = 1./cumulativeProbabilities[populationSize-1];
        for(int i=0;i<populationSize;++i){
            cumulativeProbabilities[i] *= reciprocal_sum;
        }
    }
    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::uniform_real_distribution<> dis(cumulativeProbabilities[0], cumulativeProbabilities[populationSize-1]);
    if(!winners){
        winners = new int[winnersSize];
    }
    for(int i=0;i<winnersSize;++i){
        float pick = dis(gen);
        winners[i] = std::lower_bound(cumulativeProbabilities, cumulativeProbabilities+populationSize, pick) - cumulativeProbabilities; 
    }
    delete[] cumulativeProbabilities;
}

void calculateRanks(float* fitness, int populationSize, bool maximizeFitness, int* ranksLookup){
    assert(ranksLookup);
    float *fitnessOrdered = new float[populationSize];
    std::memcpy(fitnessOrdered, fitness, populationSize*sizeof(float));
    std::sort(fitnessOrdered, fitnessOrdered+populationSize-1);
    std::unordered_map<float, int> fitnessMap;
    fitnessMap.reserve(populationSize);
    for(int i=0;i<populationSize;++i){
        fitnessMap[fitness[i]] = i;
    }
    if(maximizeFitness){
        for(int i=0;i<populationSize;++i){
            ranksLookup[populationSize - 1 - i] = fitnessMap[fitnessOrdered[i]];
        }
    } else {
        for(int i=0;i<populationSize;++i){
            ranksLookup[i] = fitnessMap[fitnessOrdered[i]];
        }
    }
    delete[] fitnessOrdered;
}

void calculateLinearRankingProbabilities(float selectionPressure, float *probabilities, int populationSize){
    assert(probabilities);
    float k2 = selectionPressure/float(populationSize - 1);
    probabilities[0] = selectionPressure;
    for(int i=1;i<populationSize-1;++i){
        probabilities[i] = probabilities[i-1] + selectionPressure - i*k2;
    }
    if(populationSize>1){
        probabilities[populationSize - 1] = probabilities[populationSize - 2];
    }
}

float* linearRankingProbabilitiesGenerator(float selectionPressure, int populationSize){
    float* temp_probabilities;
    auto temp_iterator = linearRanksData.find(selectionPressure);
    RanksData *temp_data = NULL;
    int previousRank = 0;
    if(temp_iterator!=linearRanksData.end()){
        auto temp_map = temp_iterator->second;
        auto inner_temp_iterator = temp_map.find(populationSize);
        if(inner_temp_iterator!=temp_map.end()){
            temp_data = &linearRanksData[selectionPressure][populationSize]; 
            for(const auto& m : linearRanksData){
                for(const auto& data : m.second){
                    if(data.second.usage>temp_data->usage){
                        previousRank++;
                    }
                }
            }
            temp_probabilities = linearRanksData[selectionPressure][populationSize].probabilities;
            temp_data->usage++;
        }
    }
    if(temp_data!=NULL){
        int rank=0;
        RanksData* nextData = temp_data;
        for(const auto& m : linearRanksData){
            for(const auto& data : m.second){
                if(data.second.usage>temp_data->usage){
                    rank++;
                } else {
                    if(data.second.usage>nextData->usage){
                        nextData = const_cast<RanksData*>(&(data.second));
                    }
                }
            }
        }
        if(rank<storedProbabilities){
            if(previousRank>=storedProbabilities){
                calculateLinearRankingProbabilities(selectionPressure, temp_probabilities, populationSize);
                delete[] nextData->probabilities;
            }
        }
    } else {
        calculateLinearRankingProbabilities(selectionPressure, temp_probabilities, populationSize);
        bool mustBreak = false;
        for(const auto& m : linearRanksData){
            if(mustBreak) {break;}
            for(const auto& data : m.second){
                if(data.second.usage==1){
                    auto mMap = const_cast<std::unordered_map<int, RanksData>*>(&(m.second));
                    mMap->erase(data.first);
                    RanksData mRanksData = {1, populationSize, temp_probabilities};
                    linearRanksData[selectionPressure][populationSize] = mRanksData;
                    mustBreak = true;
                    break;
                }
            }
        }
    }
    return temp_probabilities;
}

void deleteProbabilities(float* probabilities, bool linear){
    bool exists = false;
    if(linear){
        for(const auto& m : linearRanksData){
            for(const auto& data : m.second){
                if(data.second.probabilities == probabilities){
                    exists = true;
                }
            }
        }
    } else {
        for(const auto& m : exponentialRanksData){
            if(m.second.probabilities == probabilities){
                exists = true;
            }
        }
    }
    if(!exists){
        delete[] probabilities;
    }
}


void linearRanking(int populationSize, float *fitness, bool maximizeFitness, float selectionPressure, int *winners, int winnersSize){
    assert(selectionPressure>1 && selectionPressure<2 && "linearRanking: selectionPressure must be between 1 and 2 , extremes excluded.\n");
    int *ranksLookup = new int[populationSize];
    calculateRanks(fitness, populationSize, maximizeFitness, ranksLookup);
    float *cumulativeProbabilities = linearRankingProbabilitiesGenerator(selectionPressure, populationSize);
    assert(winners);
    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::uniform_real_distribution<> dis(cumulativeProbabilities[0], cumulativeProbabilities[populationSize-1]);
    for(int i=0;i<winnersSize;++i){
        winners[i] = ranksLookup[std::lower_bound(cumulativeProbabilities, cumulativeProbabilities+populationSize, dis(gen)) - cumulativeProbabilities];
    }
    deleteProbabilities(cumulativeProbabilities, true);
    delete[] ranksLookup;
}

void calculateExponentialRankingProbabilities(float k1, float *probabilities, int populationSize, float startingValue, int startingIndex){
    assert(startingIndex>=0 && "calculateExponentialRankingProbabilities: startingIndex must be positive.\n");
    assert(startingValue>=0 && "calculateExponentialRankingProbabilities: startingValue must be positive.\n");
    assert(probabilities);
    if(startingIndex) {
        probabilities[0] = startingValue + k1*std::pow(1.-k1, float(startingIndex+1));
    } else {
        probabilities[0] = k1;
    }
    for(int i=1;i<populationSize - startingIndex;i++){
        probabilities[i] = probabilities[i-1] + k1*std::pow(1.-k1, float(i+startingIndex));
    }
}

float* exponentialRankingProbabilitiesGenerator(float k1, int populationSize){
    float* temp_probabilities;
    auto temp_iterator = exponentialRanksData.find(k1);
    RanksData *temp_data = NULL;
    int previousRank = 0;
    if(temp_iterator!=exponentialRanksData.end()){
        temp_data = &(temp_iterator->second); 
            for(const auto& m : exponentialRanksData){
                if(m.second.usage>temp_data->usage){
                    previousRank++;
                }
            }
            temp_probabilities = temp_data->probabilities;
            temp_data->usage++;
        if(temp_data!=NULL){
            int rank=0;
            RanksData* nextData = temp_data;
            for(const auto& m : exponentialRanksData){
                if(m.second.usage>temp_data->usage){
                    rank++;
                } else {
                    if(m.second.usage>nextData->usage){
                        nextData = const_cast<RanksData*>(&(m.second));
                    }
                }
            }
            if(rank<storedProbabilities){
                if(previousRank>=storedProbabilities){
                    calculateExponentialRankingProbabilities(k1, temp_probabilities, populationSize, 0., 0);
                    delete[] nextData->probabilities;
                    temp_data->populationSize = populationSize;
                }
                else{
                    int currentPopSize = temp_data->populationSize;
                    if(populationSize>currentPopSize){
                        calculateExponentialRankingProbabilities(k1, temp_probabilities, populationSize, temp_probabilities[currentPopSize-1], currentPopSize-1);
                        temp_data->populationSize = populationSize;
                    }
                }
            }
        }
    } 
    if(temp_data==NULL){
        calculateExponentialRankingProbabilities(k1, temp_probabilities, populationSize, 0., 0);
        bool mustBreak = false;
        for(const auto& m : exponentialRanksData){
            if(mustBreak) {break;}
            if(m.second.usage==1){
                exponentialRanksData.erase(m.first);
                RanksData mRanksData = {1, populationSize, temp_probabilities};
                exponentialRanksData[k1] = mRanksData;
                mustBreak = true;
                break;
            }
        }
    }
    return temp_probabilities;
}


void exponentialRanking(int populationSize, float *fitness, bool maximizeFitness, float k1, int *winners, int winnersSize){
    assert(populationSize>0 && "exponentialRanking: populationSize must be positive.\n");
    assert(k1 >= 0.01 && k1 <= 0.1 && "exponentialRanking: k1 must be between 0.01 and 0.1\n");
    int *ranksLookup = new int[populationSize];
    calculateRanks(fitness, populationSize, maximizeFitness, ranksLookup);
    float* cumulativeProbabilities = exponentialRankingProbabilitiesGenerator(k1, populationSize);
    if(!winners){
        winners = new int[winnersSize];
    }
    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::uniform_real_distribution<> dis(cumulativeProbabilities[0], cumulativeProbabilities[populationSize-1]);
    for(int i=0;i<winnersSize;++i){
        winners[i] = ranksLookup[std::lower_bound(cumulativeProbabilities, cumulativeProbabilities+populationSize, dis(gen)) - cumulativeProbabilities];
    }
    deleteProbabilities(cumulativeProbabilities, false);
    delete[] ranksLookup;
}

void tournamentRanking(int populationSize, float *fitness, bool maximizeFitness, int tournamentSize, int *winners, int winnersSize){
    assert(tournamentSize > 1 && "tournamentRanking: tournamentSize must be greater than 1");
    assert(tournamentSize < populationSize && "tournamentRanking: tournamentSize must be less than populationSize");
    int *ranksLookup = new int[populationSize];
    calculateRanks(fitness, populationSize, maximizeFitness, ranksLookup);
    assert(winners);
    for(int k=0;k<winnersSize;++k){
        //Inside-out Fisher-Yattes Shuffle to get a random permutation of [0, populationSize-1]
        std::random_device rd;  
        std::mt19937 gen(rd()); 
        int permutation[populationSize];
        for(int i=0;i<populationSize;++i){
            std::uniform_int_distribution<> dis(0, i);
            int j = dis(gen);
            if(j!=i){
                permutation[i] = permutation[j];
            }
            permutation[j] = i;
        }
        winners[k] = permutation[0];
        for(int i=1;i<tournamentSize;++i){
            if(winners[k] > permutation[i]){
                winners[k] = permutation[i];
            }
        }
    }
    for(int i=0;i<winnersSize;++i){
        winners[i] = ranksLookup[winners[i]];
    }
    delete[] ranksLookup;
}

void twoPointsCrossover(uint8_t *parent1, uint8_t *parent2, int length, uint8_t *child, uint64_t *genesLoci, int genesLociLength){
   assert(length>2 && "twoPointsCrossover: can't crossover genomes of size less than 3.\n");
   assert(genesLociLength>2 && "twoPointsCrossover: can't crossover genomes with less than 3 genes.\n");
   assert(std::is_sorted(genesLoci, genesLoci+genesLociLength) && "twoPointsCrossover: genesLoci needs to be sorted in non-descending order.\n");
   assert(child);
   std::random_device rd;  
   std::mt19937 gen(rd()); 
   std::uniform_int_distribution<> dis(0, genesLociLength-1);
   uint64_t cut1 = genesLoci[dis(gen)];
   uint64_t cut2;
   do{cut2 = genesLoci[dis(gen)];} while(cut2==cut1);
   if(cut1>cut2){
        uint64_t dummy = cut2;
        cut2 = cut1;
        cut1= dummy;
   }
   std::memcpy(child, parent1, cut1);
   std::memcpy(child+cut1, parent2+cut1, cut2 - cut1);
   std::memcpy(child+cut2, parent1+cut2, length - cut2);
}

void twoPointsCrossover(uint8_t *parent1, uint8_t *parent2, int length, uint8_t *child){
   assert(length>2 && "twoPointsCrossover: can't crossover genomes of size less than 3.\n");
   assert(child);
   std::random_device rd;  
   std::mt19937 gen(rd()); 
   std::uniform_int_distribution<> dis(1, length-2);
   int cut1 = dis(gen);
   int cut2;
   do{cut2 = dis(gen);} while(cut2==cut1);
   if(cut1>cut2){
        int dummy = cut2;
        cut2 = cut1;
        cut1= dummy;
   }
   std::memcpy(child, parent1, cut1);
   std::memcpy(child+cut1, parent2+cut1, cut2 - cut1);
   std::memcpy(child+cut2, parent1+cut2, length - cut2);
}


void uniformCrossover(uint8_t *parent1, uint8_t *parent2, int length, uint8_t *child, uint64_t *genesLoci, int genesLociLength){
   assert(std::is_sorted(genesLoci, genesLoci+genesLociLength) && "uniformCrossover: genesLoci needs to be sorted in non-descending order.\n");
    assert(child);
    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::uniform_int_distribution<uint64_t> dis(0, UINT64_MAX);
    uint64_t mask = dis(gen);
    if(genesLoci[0]!=0){
        uint64_t *_genesLoci = new uint64_t[genesLociLength+2];
        _genesLoci[0] = 0;
        _genesLoci[genesLociLength+1] = length;
        std::memcpy(_genesLoci+1, genesLoci, genesLociLength);
        for(int i=0;i<genesLociLength+2;++i){
            uint8_t maskBit = mask & (1<<i);
            for(uint64_t j=_genesLoci[i];j<_genesLoci[i+1];++j){
                child[i] = parent1[i]*maskBit + parent2[i]*(1 - maskBit);
            }
        }
        delete[] _genesLoci;
    } else {
        uint64_t *_genesLoci = new uint64_t[genesLociLength+1];
        _genesLoci[genesLociLength] = length;
        std::memcpy(_genesLoci, genesLoci, genesLociLength);
        for(int i=0;i<genesLociLength+1;++i){
            uint8_t maskBit = mask & (1<<i);
            for(uint64_t j=_genesLoci[i];j<_genesLoci[i+1];++j){
                child[i] = parent1[i]*maskBit + parent2[i]*(1 - maskBit);
            }
        }
        delete[] _genesLoci;
    }
}

void uniformCrossover(uint8_t *parent1, uint8_t *parent2, int length, uint8_t *child){
    assert(child);
    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::uniform_int_distribution<uint64_t> dis(0, UINT64_MAX);
    uint64_t mask = dis(gen);
    for(int i=0;i<length;++i){
        uint8_t maskBit = mask & (1<<i);
        child[i] = parent1[i]*maskBit + parent2[i]*(1 - maskBit);
    }
}

void mutate(uint8_t *individual, int length, float mutationProbability){
    assert(mutationProbability>0 && "mutate: mutationProbability must be greater than 0.\n");
    assert(mutationProbability<1 && "mutate: mutationProbability must be less than 1.\n");
    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::uniform_real_distribution<float> prob(0, 1);
    std::uniform_int_distribution<uint8_t> dis8(0, 7);
    for(int i=0;i<length;++i){
        if(prob(gen)<mutationProbability){
            *individual^= (1<<dis8(gen));
        }
    }
}
