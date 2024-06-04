#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <map>
#include <utility>
#include <iostream>

namespace Utils {

    std::pair<int, std::map<int, int>> CountUniqueElements(const std::vector<int>& vec) {
        std::map<int, int> frequencyMap;
        for (int element : vec) {
            frequencyMap[element]++;
        }
        int uniqueCount = frequencyMap.size();
        return std::make_pair(uniqueCount, frequencyMap);
    }

    int PrintElementsWithMinFrequency(const std::map<int, int>& frequencyMap, int minFrequency) {
        int count = 0;
        for (const auto& pair : frequencyMap) {
            if (pair.second >= minFrequency) {
                count++;
            }
        }

        return count;
        //std::cout << "Number of unique elements with frequency >= " << minFrequency << ": " << count << std::endl;

        //std::cout << "Elements with frequency >= " << minFrequency << ":" << std::endl;
        //for (const auto& pair : frequencyMap) {
        //    if (pair.second >= minFrequency) {
        //        std::cout << pair.first << ": " << pair.second << std::endl;
        //    }
        //}
    }

    std::vector<int> DivideElementsBy(std::vector<int>& vec, int divisor) {
        std::vector<int> result;
        result.reserve(vec.size());
        for (int element : vec) {
            result.push_back(element / divisor);
        }
        return result;
    }

}

#endif // UTILS_HPP

