#pragma once

#include <vector>
#include <stdexcept>

std::vector<int> Unique(const std::vector<int>& data) {
    std::vector<int> ans;
    if (data.empty()) {
        return ans;
    }
    ans.push_back(data[0]);
    for (std::size_t i = 1; i < data.size(); ++i) {
        if (data[i] != data[i - 1]) {
            ans.push_back(data[i]);
        }
    }
    return ans;
}
