#ifndef NAMESMAP_H
#define NAMESMAP_H
#pragma once
#include <unordered_map>
#include <vector>
#include <string>
/**
 * @brief A class that maps names to indices and vice versa.
 *
 * The NamesMap class provides a bidirectional mapping between strings (names)
 * and integers (indices). It can be used to map row or column names to indices
 * in a matrix or data structure.
 */
class NamesMap
{
public:
    NamesMap() {}
    /**
     * @brief Adds a name and associates it with an index.
     *
     * This method stores the mapping between the provided name and index
     * in both directions (name to index and index to name).
     *
     * @param name The name to be associated with the index.
     * @param index The index to associate with the name.
     */
    void add(const std::string &name, int index);
    /**
     * @brief Gets the index associated with a name.
     *
     * This method returns the index associated with the provided name.
     * If the name is not found, it returns -1.
     *
     * @param name The name to lookup.
     * @return The index associated with the name, or -1 if not found.
     */
    int get_index(const std::string &name) const;
    /**
     * @brief Gets the name associated with an index.
     *
     * This method returns the name associated with the provided index.
     * If the index is not found, it returns an empty string.
     *
     * @param index The index to lookup.
     * @return The name associated with the index, or an empty string if not found.
     */
    std::string get_name(int index) const;
    /**
     * @brief Gets a vector of all names in the map.
     *
     * This method returns a vector containing all names currently stored
     * in the map, in no specific order.
     *
     * @return A vector of strings containing all names.
     */
    std::vector<std::string> get_all_names() const;

private:
    std::unordered_map<std::string, int> name_to_index;
    std::unordered_map<int, std::string> index_to_name;
};
#endif