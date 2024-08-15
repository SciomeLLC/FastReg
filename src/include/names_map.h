#ifndef NAMESMAP_H
#define NAMESMAP_H
#pragma once
#include <unordered_map>
#include <vector>
#include <string>

// Custom data structure to store row or column names
class NamesMap {
public:
    NamesMap() {}
    void add(const std::string& name, int index);
    int get_index(const std::string& name) const;
    std::string get_name(int index) const;
    std::vector<std::string> get_all_names() const;

private:
    std::unordered_map<std::string, int> name_to_index;
    std::unordered_map<int, std::string> index_to_name;
};
#endif