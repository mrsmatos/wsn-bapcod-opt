/**
 *
 * This file bcIndexC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef INDEXC_H
#define INDEXC_H

#include "bcUsefulHeadFil.hpp"
#include "bcMultiIndexC.hpp"

#if(BOOST_VERSION >= 104000)
#include <boost/unordered_map.hpp>
#else
#include <unordered_map>
#endif //BOOST_LIB_VERSION


struct Index
{
  int _ref;

  Index(): _ref(0) {}

  Index(const Index & that):  _ref(that._ref) {}
  
  Index(const int & ref): _ref(ref) {}

  virtual ~Index() {}

  virtual const int & ref() const
  {
    return _ref;
  }

  virtual void ref(const int & r) 
  {
    _ref = r;
  }

  /// Operator conversion:
  operator int() const  
  {
    return _ref;
  }
  
  const Index operator-() const 
  {
    return -_ref;
  }

  virtual std::ostream & print(std::ostream& os = std::cout) const
  {
    os   << "ref = " << _ref << std::endl;
	return os;
  }
};



template <typename ObjectType>
class MapByMultiIndex: public std::map < MultiIndex , ObjectType >
{
};

struct ihash
    : std::unary_function<MultiIndex, std::size_t>
{
    std::size_t operator()(MultiIndex const& e) const;
};

struct iequal_to
    : std::binary_function<MultiIndex, MultiIndex, bool>
{
    bool operator()(MultiIndex const& x, MultiIndex const& y) const
    {
        return (x==y);
    }
};


#if (BOOST_VERSION >= 104000)
template <typename ObjectType>
class UnorderedMapByMultiIndex: public boost::unordered_map< MultiIndex , ObjectType, ihash>
{
};
#else
#ifndef _MSC_VER
// Only available on Unix
template <typename ObjectType>
class UnorderedMapByMultiIndex: public std::unordered_map< MultiIndex , ObjectType, ihash>
{
};
#else
// Only available on Windows
template <typename ObjectType>
class UnorderedMapByMultiIndex: public stdext::hash_map< MultiIndex , ObjectType>
{
};
#endif //_MSC_VER
#endif

template <typename ObjectType>
class MapByNameAndMultiIndex: public std::map < std::string, std::map < MultiIndex , ObjectType > >
{
};


class ColGenSpConf;
typedef MapByMultiIndex< ColGenSpConf * > MapColGenSpConfPtrByMultiIndex;

class ColGenSpConf;
typedef MapByNameAndMultiIndex< ColGenSpConf * > MapColGenSpConfPtrByNameAndMultiIndex;

class ProbConfig;

class InstanciatedVar;
class NetworkArc;
class NetworkVertex;
typedef MapByMultiIndex< InstanciatedVar * > MapInstVarPtrByMultiIndex;

class InstanciatedConstr;

typedef UnorderedMapByMultiIndex< InstanciatedConstr * > MapInstConstrPtrByMultiIndex;

class Double;
typedef MapByMultiIndex< Double > MapValByMultiIndex;


template <typename PointerType>
struct SmartIndex: public Index
{
  PointerType * _propertyClassPtr;

 SmartIndex(): Index(), _propertyClassPtr(NULL) {}

 SmartIndex(PointerType * propertyClassPtr): Index(), _propertyClassPtr(propertyClassPtr) {}

  PointerType * ptr() {return _propertyClassPtr;}

};

/// Comparators <
inline bool operator<(const Index & x, const Index & y) 
{
  return (x._ref < y._ref);
}

inline bool operator<(const Index & x, const int & y) 
{
  return (x._ref < y);
}

inline bool operator<(const int & x, const Index & y) 
{
  return (x < y._ref);
}

/// Comparators <=
inline bool operator<=(const Index & x, const Index & y) 
{
  return (x._ref <= y._ref);
}
inline bool operator<=(const Index & x, const int & y) 
{
  return (x._ref <= y);
}

inline bool operator<=(const int & x, const Index & y) 
{
  return (x <= y._ref);
}

/// Comparators >
inline bool operator>(const Index & x, const Index & y) 
{
  return (x._ref > y._ref);
}

inline bool operator>(const Index & x, const int & y) 
{
  return (x._ref > y);
}

inline bool operator>(const int & x, const Index & y) 
{
  return (x > y._ref);
}

/// Comparators >=
inline bool operator>=(const Index & x, const Index & y) 
{
  return (x._ref >= y._ref);
}

inline bool operator>=(const Index & x, const int & y) 
{
  return (x._ref >= y);
}

inline bool operator>=(const int & x, const Index & y) 
{
  return (x >= y._ref);
}

/// Comparators ==
inline bool operator==(const Index & x, const Index & y) 
{
  return (x._ref == y._ref);
}

inline bool operator==(const Index & x, const int & y) 
{
  return (x._ref == y);
}

inline bool operator==(const int & x, const Index & y) 
{
  return (x == y._ref);
}

/// Comparators !=
inline bool operator!=(const Index & x, const Index & y) 
{
  return (x._ref != y._ref);
}

inline bool operator!=(const Index & x, const int & y) 
{
  return (x._ref != y);
}

inline bool operator!=(const int & x, const Index & y) 
{
  return (x != y._ref);
}

/// Operator << >>
inline std::ostream& operator<<(std::ostream& os, const Index & x) 
{
  return os << x._ref;
}

inline std::istream& operator>>(std::istream& is, Index& x) 
{
  is >>  x._ref;
  return is; 
}

/// Other operators int
inline const Index operator+(const Index & x, const int & y) 
{
  return (x._ref +  y);
}

inline const Index operator+(const int & x, const Index & y) 
{
  return (x + y._ref);
}

inline const Index operator+(const Index & x, const Index & y) 
{
  return (x._ref + y._ref);
}

inline const Index operator-(const Index & x, const int & y) 
{
  return (x._ref - y);
}

inline const Index operator-(const int & x, const Index & y) 
{
  return (x - y._ref);
}

inline const Index operator-(const Index & x, const Index & y) 
{
  return (x._ref - y._ref);
}

inline const Index operator*(const Index & x, const int & y) 
{
  return (x._ref *  y);
}

inline const Index operator*(const int & x, const Index & y) 
{
  return (x * y._ref);
}

inline const Index operator*(const Index & x, const Index & y) 
{
  return (x._ref * y._ref);
}

inline const Index operator/(const Index & x, const int & y) 
{
  return (x._ref / y);
}

inline const Index operator/(const int & x, const Index & y) 
{
  return ( x / y._ref);
}

inline const Index operator/(const Index & x, const Index & y) 
{
  return ( x._ref / y._ref);
}

inline Index & operator+=(Index & x, const int & y) 
{
  ((x._ref) +=  y); 
  return x;
}

inline int & operator+=(int & x, const Index & y) 
{
  return (x +=  y._ref);
}

inline Index & operator+=(Index & x, const Index & y) 
{
  (x._ref +=  y._ref); 
  return x;
}

inline Index & operator-=(Index & x, const int & y) 
{
  (x._ref -=  y); 
  return x;
}

inline int & operator-=(int & x, const Index & y) 
{
  return ( x -= y._ref);
}

inline Index & operator-=(Index & x, const Index & y) 
{
  (x._ref -= y._ref); 
  return x;
}

inline Index & operator*=(Index & x, const int & y) 
{
  (x._ref *=  y); 
  return x;
}

inline int & operator*=(int & x, const Index & y) 
{
  return (x *= y._ref);
}

inline Index & operator*=(Index & x, const Index & y) 
{
  (x._ref *= y._ref); 
  return x;
}

inline Index & operator/=(Index & x, const int & y) 
{
  x._ref /=  y; 
  return x;
}

inline int & operator/=(int & x, const Index & y) 
{
  return (x =  x / y._ref);
}

inline Index & operator/=(Index & x, const Index & y) 
{
  (x._ref /= y._ref); 
  return x;
}

struct IndexPtrSort
{
  bool operator()(Index * x, Index * y) const 
  {
    return (x->_ref < y->_ref);
  }

};


template < typename PointerType >
class SmartIndexRange: public std::vector< SmartIndex < PointerType > >
{
  int _lowerBound;
  int _upperBound;
  int _last;

 public: 
 SmartIndexRange():
  _lowerBound(0),
    _upperBound(0),
    _last(0)
      {}

 SmartIndexRange(const std::vector< PointerType > & vect) :
  std::vector< SmartIndex < PointerType > >(vect.size()),
    _lowerBound(0),
    _upperBound(vect.size()),
    _last(vect.size()-1)
      {
	for (typename std::vector< PointerType >::iterator it = vect.begin(); 
	     it != vect.end(); ++it) this->push_back(SmartIndex < PointerType > (*it));
      }

  void add (PointerType * pointer) 
  {
    this->push_back(SmartIndex < PointerType > (pointer));
    ++_upperBound;
    _last = _upperBound - 1;
    return;
  }

  SmartIndex < PointerType > & bringIntoRange(const int & ref)
    {
      if (ref < _lowerBound) return (*this)[0];
      if (ref > _upperBound) return (*this)[_last];

      return (*this)[ref - _lowerBound];
    }

  const bool check(const int & ref)
  {
    if (ref < _lowerBound) return true;
    if (ref > _upperBound) return true;
    return false;
  }

  SmartIndex < PointerType > & find(const int & ref)
    {
      if (check(ref)) 
	{
	  std::cout << " Index out of Range " << std::endl;
	  exit(1);
	}

      return (*this)[ref - _lowerBound];
    }


};

template < typename FirstPointerType, typename SecondPointerType, typename ThirdPointerType, typename FourthPointerType, typename FifthPointerType, typename SixthPointerType >  
class SmartIndexCell : public MultiIndex
{
  bool _simpleMutliIndex;
  FirstPointerType * _firstPtr;
  SecondPointerType * _secondPtr;
  ThirdPointerType * _thirdPtr;
  FourthPointerType * _fourthPtr;
  FifthPointerType * _fifthPtr;
  SixthPointerType * _sixthPtr;
  
 public:

 SmartIndexCell(const MultiIndex & multiIndex = MultiIndex()):
  MultiIndex(multiIndex),
    _simpleMutliIndex(true),
    _firstPtr(NULL),
    _secondPtr(NULL),
    _thirdPtr(NULL),
    _fourthPtr(NULL),
    _fifthPtr(NULL),
    _sixthPtr(NULL)
      {
      }

 SmartIndexCell(FirstPointerType * i1, 
		SecondPointerType * i2 = NULL, 
		ThirdPointerType * i3 = NULL, 
		FourthPointerType * i4 = NULL, 
		FifthPointerType * i5 = NULL, 
		SixthPointerType * i6 = NULL):
  MultiIndex((i1 != NULL? i1->ref() : -1),
	     (i2 != NULL? i2->ref() : -1),
	     (i3 != NULL? i3->ref() : -1),
	     (i4 != NULL? i4->ref() : -1),
	     (i5 != NULL? i5->ref() : -1),
	     (i6 != NULL? i6->ref() : -1)),
    _simpleMutliIndex(false),
    _firstPtr(i1), 
    _secondPtr(i2), 
    _thirdPtr(i3), 
    _fourthPtr(i4), 
    _fifthPtr(i5), 
    _sixthPtr(i6) 
    {
    }
  
  FirstPointerType * firstPtr() const {return _firstPtr;}
  SecondPointerType * secondPtr() const {return _secondPtr;}
  ThirdPointerType * thirdPtr() const {return _thirdPtr;}
  FourthPointerType * fourthPtr() const {return _fourthPtr;}
  FifthPointerType * fifthPtr() const {return _fifthPtr;}
  SixthPointerType * sixthPtr() const {return _sixthPtr;}



  virtual ~SmartIndexCell()
    {
    }
};



/// DEFAUT TEMPLATE ARGUMENTS

template < typename FirstPointerType, typename SecondPointerType, typename ThirdPointerType, typename FourthPointerType, typename FifthPointerType, typename SixthPointerType = Index > class SmartIndexCell;

template < typename FirstPointerType, typename SecondPointerType, typename ThirdPointerType, typename FourthPointerType, typename FifthPointerType = Index, typename SixthPointerType > class SmartIndexCell;

template < typename FirstPointerType, typename SecondPointerType, typename ThirdPointerType, typename FourthPointerType = Index, typename FifthPointerType, typename SixthPointerType > class SmartIndexCell;

template < typename FirstPointerType, typename SecondPointerType, typename ThirdPointerType = Index, typename FourthPointerType, typename FifthPointerType, typename SixthPointerType > class SmartIndexCell;

template < typename FirstPointerType, typename SecondPointerType = Index, typename ThirdPointerType, typename FourthPointerType, typename FifthPointerType, typename SixthPointerType > class SmartIndexCell;

template < typename FirstPointerType = Index, typename SecondPointerType, typename ThirdPointerType, typename FourthPointerType, typename FifthPointerType, typename SixthPointerType > class SmartIndexCell;


typedef SmartIndexCell<Index,Index,Index,Index,Index,Index> IndexCell;


struct IndexCellSort
{
  bool operator()(const IndexCell & a, const IndexCell & b) const
  {
    return (a < b);
  }
};


template <typename PointerType>
class DirectAccesMapByIndex
{
  PointerType ** _1indexArray;
  PointerType *** _2indexArray;
  PointerType **** _3indexArray;
  PointerType ***** _4indexArray;
  PointerType ****** _5indexArray;
  PointerType ******* _6indexArray;
  int _sizeFirst;
  int _sizeSecond;
  int _sizeThird;
  int _sizeFourth;
  int _sizeFith;
  int _sizeSixth;

 public: 
 DirectAccesMapByIndex():
  _1indexArray(NULL), 
    _2indexArray(NULL),
    _3indexArray(NULL),
    _4indexArray(NULL),
    _5indexArray(NULL),
    _6indexArray(NULL),
    _sizeFirst(0),
    _sizeSecond(0),
    _sizeThird(0),
    _sizeFourth(0),
    _sizeFith(0),
    _sizeSixth(0)
      {}
  void initSize(const int & sizeFirst)
    {
      _sizeFirst = sizeFirst;
      _1indexArray = new PointerType*[_sizeFirst];
      for (int i = 0; i < _sizeFirst; ++i) _1indexArray[i] = NULL;
      return;
    }

  void initSize(const int & sizeFirst, const int & sizeSecond)
    {
      _sizeFirst = sizeFirst;
      _sizeSecond = sizeSecond;
      _2indexArray = new PointerType**[_sizeFirst];
      for (int i = 0; i < _sizeFirst; ++i) 
	{
	  _2indexArray[i] = new PointerType*[_sizeSecond];
	  for (int j = 0; j < _sizeSecond; ++j) 
	    {
	      _2indexArray[i][j] = NULL;
	    }
	}
      return;
    }

  PointerType * find(const Index & first)
  {
    if (_1indexArray == NULL) return NULL;
    if (first >= _sizeFirst) exit(1);
    return _1indexArray[first];
  }

  PointerType * find(const Index & first, const Index & second)
  {
    if (_2indexArray == NULL) return NULL;
    if (first >= _sizeFirst) exit(1);
    if (second >= _sizeSecond) exit(1);
    return _2indexArray[first][second];
  }

  ~DirectAccesMapByIndex()
    {
      if (_1indexArray != NULL) delete _1indexArray;
      if (_2indexArray != NULL) 
	for (int i = 0; i < _sizeFirst; ++i)  delete _2indexArray[i];

    }
};



#endif // INDEXC_H
