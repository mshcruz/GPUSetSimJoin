#include <algorithm> 
#include <cctype>    
#include <chrono>
#include <cstdlib>
#include <functional>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <regex>
#include <sstream>
#include <string>
#include <sys/time.h>
#include <unordered_map>
#include <vector>

typedef unsigned int uint;
typedef unsigned long ulong;

int fastMod(uint a) { //Mod 65535
  a = (a >> 16) + (a & 0xFFFF); /* sum base 2**16 digits */
  if (a < 65535) return a;
  if (a < (2 * 65535)) return a - 65535;
  return a - (2 * 65535);
}

template< typename K, typename T, typename C, typename A >
  typename std::multimap<K,T,C,A>::size_type num_unique_keys( const std::multimap<K,T,C,A>& mmap )
{
  if( mmap.size() < 2U ) return mmap.size() ;
  const C& cmp = mmap.key_comp() ;
  typename std::multimap<K,T,C,A>::size_type n = 1 ;
  auto prev = mmap.begin()->first ;
  for( auto iter = mmap.begin() ; iter != mmap.end() ; ++iter )
    if( cmp( prev, iter->first ) ) { ++n ; prev = iter->first ; }
  return n ;
}

std::string 
removeNonAlphaNum(std::string str)
{
  for (uint i = 0; i < str.size(); i++) {
    if (!isalnum(str[i])) {
      if (str[i] != ' ') {
	str.erase(i,1);	
      }
    }
  }
  return str;
}

class Timer
{
public:
    void start() { m_start = my_clock(); }
    void stop() { m_stop = my_clock(); }
    double elapsed_time() const {
        return m_stop - m_start;
    }

private:
    double m_start, m_stop;
    double my_clock() const {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        return tv.tv_sec + tv.tv_usec * 1e-6;
    }
};

std::vector<int>
decimalToBinaryPositions(int n)
{
  std::vector<int> binaryPositions;
  for (int i = 0; n > 0; n = n >> 1, i++) {
    if (n%2) {
      binaryPositions.push_back(i);
    }
  }
  return binaryPositions;
}

std::unordered_map<std::string, int>
createStopWordsMap()
{
  std::vector<std::string> stopWords = {"a","able","about","across","after","all","almost","also","am","among","an","and","any","are","as","at","be","because","been","but","by","can","cannot","could","dear","did","do","does","either","else","ever","every","for","from","get","got","had","has","have","he","her","hers","him","his","how","however","i","if","in","into","is","it","its","just","least","let","like","likely","may","me","might","most","must","my","neither","no","nor","not","of","off","often","on","only","or","other","our","own","rather","said","say","says","she","should","since","so","some","than","that","the","their","them","then","there","these","they","this","tis","to","too","twas","us","wants","was","we","were","what","when","where","which","while","who","whom","why","will","with","would","yet","you","your"};
  std::unordered_map<std::string, int> stopWordsMap;
  
  for (uint i = 0; i < stopWords.size(); i++) {
    std::pair<std::string, int> stopWord(stopWords[i],i);
    stopWordsMap.insert(stopWord);
  }

  return stopWordsMap;
}

std::string
stringToLower(std::string str)
{
  for (uint i = 0; i < str.size(); i++) {
    str[i] = tolower(str[i]);
  }
  return str;
}

int //Modified from boost library
myHashing(int valueToBeHashed)
{
  int seed = 0;
  for(int i = 0; i < 7; i++)
    {
      seed ^= (int) (valueToBeHashed >> i) + (seed<<6) + (seed>>2);
    }
  seed ^= (int) valueToBeHashed + (seed<<6) + (seed>>2);
  return abs(seed);
}

inline void
hash_combine(int& seed, int value)
{
    seed ^= (int)abs(myHashing(value) + 0x9e3779b9 + (seed<<6) + (seed>>2));
}

inline int //Modified from boost library
myUniversalHashing(int a, int b, int prime, int limit, int valueToBeHashed)
{
  return ((((a*valueToBeHashed)+b)%prime)%limit);
}

int
processInputRelation(std::multimap<int,int>& shingleSetMap, std::map<std::string,int>& shingles, std::vector<std::string>& setsIDs, std::string fileName, int relationOffset)
{
  std::ifstream relation (fileName);
  std::string line = "";
  std::unordered_map<std::string,int> stopWordsMap = createStopWordsMap();
  int shingle, relationSize = 0;
  std::map<std::string,int>::iterator shingleIterator;
  if (relation.is_open()){
    while (getline(relation,line)){
      std::string value = "";
      std::string key = "";
      std::string word = "";
      std::istringstream tuple(line); 

      //Obtain the ID (key) of the record an add to the list of IDs
      getline(tuple,key,'\t');


      //Obtain the value of the record and clean the string
      getline(tuple,value,'\t');
      if (value.length() < 1) {
	std::cout << "Empty Record Error" << std::endl;
	exit(-1);
      }

      value = removeNonAlphaNum(value);
      /* std::cout << "value before: " << value << std::endl; */ //Only works for g++ 4.9 or later
      /* value = std::regex_replace(value, std::regex(" +"), std::string(" ")); */

      //Remove stop words
      std::istringstream sentence(value);
      value = "";
      while(getline(sentence,word,' ')){
	word = stringToLower(word);
	if (stopWordsMap.find(word) == stopWordsMap.end() && word.length() > 1) {
	  value += word + " ";
	}
      }

      //If the string is valid after the stop words removal, each new word becomes a key of the ordered map and it has one or more sets associated to it
      if (value.length() > 1) {
	std::istringstream sentence(value);
	while(getline(sentence,word,' ')){
	  word = stringToLower(word);
	  shingleIterator = shingles.find(word);
	  if (shingleIterator == shingles.end()) {
	    shingle = shingles.size();
	    std::pair<std::string,int> shingleIdxPair(word,shingles.size());
	    shingles.insert(shingleIdxPair);
	  } else {
	    shingle = shingleIterator -> second;
	  }
	  std::pair<int,int> shingleSetPair (relationOffset,shingle);
	  //	  std::cout << "shingleSetPair - set: " << relationOffset << " | shingle: " << shingle << " | word: " << word << "\n";
	  shingleSetMap.insert(shingleSetPair);
	}
	setsIDs.push_back(key);
	relationSize++;
	relationOffset++;
      }
    }    
    relation.close();
  } else {
    std::cout << "Error opening file.";
  }

  return relationSize;
}

template<typename IntType>
class Range {
public:
    constexpr Range(IntType e): _start(0), _end(e), _step(1) { }
    constexpr Range(IntType s, IntType e): _start(std::min(s, e)), _end(e), _step(1) { }
    // It's invalid to pass arguments such that
    // IntType is unsigned and e < s (i.e., e - s becomes negative)
    constexpr Range(IntType s, IntType e, IntType p):
        _start(s), _end((e - s) / p > 0 ? ((e - s) / p + ((e - s) % p != 0)) * p + s : s), _step(p) { }

    class iterator: public std::iterator<std::input_iterator_tag, IntType> {
    public:
        constexpr iterator(IntType s, const Range &parent): cur(s), parent(parent) { }
        constexpr IntType operator*() const { return cur; }
        iterator &operator++() { cur += parent._step; return *this; };
        const iterator operator++(int) {
            cur += parent._step;
            return iterator(cur - parent._step, parent);
        }

        constexpr bool operator==(const iterator &rhs) { return cur == rhs.cur; }
        constexpr bool operator!=(const iterator &rhs) { return cur != rhs.cur; }

    private:
        IntType cur;
        const Range &parent;
    };

    constexpr iterator begin() { return iterator(_start, *this); }
    constexpr iterator end() { return iterator(_end, *this); }

private:
    const IntType _start;
    const IntType _end;
    const IntType _step;
};

// Convinience functions for Range
template<typename IntType>
constexpr inline Range<IntType> range(IntType e)
{
    return Range<IntType>(e);
}

template<typename IntType>
constexpr inline Range<IntType> range(IntType s, IntType e)
{
    return Range<IntType>(s, e);
}

template<typename IntType>
constexpr inline Range<IntType> range(IntType s, IntType e, IntType p)
{
    return Range<IntType>(s, e, p);
}
