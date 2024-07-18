/**
 *
 * This file bcVcIndexInProbC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCVARCONSTRINDEXPROBLEMC_HPP_
#define BCVARCONSTRINDEXPROBLEMC_HPP_

#include <vector>
#include <set>
#include <sstream>

#include "bcVarConstrC.hpp"
#include "bcVcIdentifierC.hpp"
#include "bcGlobalException.hpp"

struct LexicographicSorting;

/// hollow vector allowp
///  insert in constant time
///  remove in constant time
///  begin in constant time
///  end in constant time
///  next in constant time
///  find in linear time

/// sorted hollow vector allow
///  insert in log time
///  remove in log time
///  begin in constant time
///  end in constant time
///  next in constant time
///  find in log time

/// full dimentional vector allow
///  insert in constant time
///  remove in  constant time
///  begin in constant time
///  end in constant time
///  next in constant time
///  find by index in constant time
///  find by nature in linear time

/// membership containers are stored as hollow vector

/// genVarConstr pool containers are stored in full dimentional vector

/// MasterCol pool containers are stored sorted hollow vector

/// MasterCut pool containers are stored sorted hollow vector
namespace VcIndexStatus
{

  /**
   * Enumaration of all possible status.
   */
  enum VcStatus
  {
    Undefined = -1, /// this id is available for re-assignement
    Active, /// currently considered as in the active problem (i.e. InProb), but not necessarily transfered to the active formulation (i.e. inForm can be true or false)
    Inactive, /// valid for the current problem, but currently set aside
    Unsuitable, /// not valid for the current problem, but currently set aside
    Infeasible, /// the bounds on this VC makes the problem infeasible
    ToBeActivated, /// temporaty status to indciate that a VarConstr needs to be activated
    Incompatible /// will be activated only if generated.       
  } ;
}

/**
 * Id associated to a VarConstr;
 */
class VcIndexInProb
{
protected:
  long _vcIndexInProb;
  VcIndexStatus::VcStatus _vcIndexStatus;

public:

  VcIndexInProb() :
  _vcIndexInProb(-1), _vcIndexStatus(VcIndexStatus::Undefined)
  {
  }

  const VcIndexStatus::VcStatus& vcIndexStatus() const
  {
    return _vcIndexStatus;
  }

  void vcIndexStatus(const VcIndexStatus::VcStatus& vcIndexStatus)
  {
    this->_vcIndexStatus = vcIndexStatus;
  }

  bool operator==(const VcIndexInProb & varR) const
  {
    return _vcIndexInProb == varR._vcIndexInProb;
  }

  bool operator!=(const VcIndexInProb & varR) const
  {
    return _vcIndexInProb != varR._vcIndexInProb;
  }

  bool operator<(const VcIndexInProb & varR) const
  {
    return _vcIndexInProb < varR._vcIndexInProb;
  }

  const long& vcIndexInProb() const
  {
    return _vcIndexInProb;
  }

  void vcIndexInProb(long id)
  {
    _vcIndexInProb = id;
  }
} ;

/**
 * Describes the neighbored for a VarConstr: VCTYPE is inherited from VcIndexInProb
 */
template<typename VCTYPE>
class VarConstrIndexPlaceHolder
{
private:
  VCTYPE * _vcPtr;

  VarConstrIndexPlaceHolder<VCTYPE> * _vcPrev;
  VarConstrIndexPlaceHolder<VCTYPE> * _vcNext;

  long _vcIndexInProb;

public:

  VarConstrIndexPlaceHolder(VCTYPE* vcPtr, long vcIndexInProb,
                            VarConstrIndexPlaceHolder<VCTYPE> * vcPrev = NULL,
                            VarConstrIndexPlaceHolder<VCTYPE> * vcNext = NULL) :
  _vcPtr(vcPtr), _vcPrev(vcPrev), _vcNext(vcNext), _vcIndexInProb(vcIndexInProb)
  {
  }

  VarConstrIndexPlaceHolder(const VarConstrIndexPlaceHolder<VCTYPE> & that) :
  _vcPtr(that._vcPtr), _vcPrev(that._vcPrev), _vcNext(that._vcNext),
  _vcIndexInProb(that._vcIndexInProb)
  {
  }

  virtual ~VarConstrIndexPlaceHolder()
  {
    _vcPrev = NULL;
    _vcNext = NULL;
    _vcPtr = NULL;
  }

  VCTYPE * vcPtr() const
  {
    return _vcPtr;
  }

  void vcPtr(VCTYPE * vcPtr)
  {
    _vcPtr = vcPtr;
  }

  VarConstrIndexPlaceHolder<VCTYPE> * nextVcInProb() const
  {
    return _vcNext;
  }

  void nextVcInProb(VarConstrIndexPlaceHolder<VCTYPE> * nextVcInProb)
  {
    this->_vcNext = nextVcInProb;
  }

  VarConstrIndexPlaceHolder<VCTYPE> * prevVcInProb() const
  {
    return _vcPrev;
  }

  void prevVcInProb(VarConstrIndexPlaceHolder<VCTYPE> * prevVcInProb)
  {
    this->_vcPrev = prevVcInProb;
  }

  long vcIndexInProb() const
  {
    return _vcIndexInProb;
  }

  void vcIndexInProb(long vcIndexInProb)
  {
    _vcIndexInProb = vcIndexInProb;
  }
} ;

/**
 * Simple stucture to chose according to the flag the type IsTrue or IsFalse.
 */
template<bool flag, class IsTrue, class IsFalse>
struct choose;

template<class IsTrue, class IsFalse>
struct choose<true, IsTrue, IsFalse>
{
  typedef IsTrue type;
} ;

template<class IsTrue, class IsFalse>
struct choose<false, IsTrue, IsFalse>
{
  typedef IsFalse type;
} ;

/**
 * Iterator template for VarConstrPlaceHolder.
 */
template<class T, bool isconst = false >
struct VarConstrIndexPlaceHolder_iterator
{
  typedef T value_type;
  typedef std::ptrdiff_t difference_type;
  typedef typename choose<isconst, const T&, T&>::type reference;
  typedef T* pointer;

  typedef typename choose<isconst, const VarConstrIndexPlaceHolder<T>*,
  VarConstrIndexPlaceHolder<T>*>::type nodeptr;

  nodeptr _vcPhPtr;

  VarConstrIndexPlaceHolder_iterator(nodeptr vcPhPtr = NULL) :
    _vcPhPtr(vcPhPtr)
  {
  }

  VarConstrIndexPlaceHolder_iterator(typename std::vector<VarConstrIndexPlaceHolder<T> *>::iterator it) :
    _vcPhPtr(*it)
  {
  }

  VarConstrIndexPlaceHolder_iterator(typename std::vector<VarConstrIndexPlaceHolder<T> *>::const_iterator it) :
    _vcPhPtr(*it)
  {
  }

  VarConstrIndexPlaceHolder_iterator(const VarConstrIndexPlaceHolder_iterator<T, false>& it) :
    _vcPhPtr(it._vcPhPtr)
  {
  }

  pointer operator*() const
  {
    return _vcPhPtr->vcPtr();
  }

  pointer operator->() const
  {
    return _vcPhPtr->vcPtr();
  }

  VarConstrIndexPlaceHolder_iterator& operator++()
  {
    _vcPhPtr = _vcPhPtr->nextVcInProb();
    return *this;
  }

  VarConstrIndexPlaceHolder_iterator operator++(int)
  {
    VarConstrIndexPlaceHolder_iterator tmp(*this);
    ++ * this;
    return tmp;
  }

  VarConstrIndexPlaceHolder_iterator& operator--()
  {
    _vcPhPtr = _vcPhPtr->prevVcInProb();
    return *this;
  }

  VarConstrIndexPlaceHolder_iterator operator--(int)
  {
    VarConstrIndexPlaceHolder_iterator tmp(*this);
    -- * this;
    return tmp;
  }

  VarConstrIndexPlaceHolder_iterator operator+=(int i)
  {
    while (i > 0)
    {
      ++ * this;
      i--;
    }
    return *this;
  }

  VarConstrIndexPlaceHolder_iterator operator+(int i)
  {
    VarConstrIndexPlaceHolder_iterator tmp(*this);
    while (i > 0)
    {
      ++tmp;
      i--;
    }
    return tmp;
  }

  VarConstrIndexPlaceHolder_iterator operator-=(int i)
  {
    while (i > 0)
    {
      -- * this;
      i--;
    }
    return *this;
  }

  VarConstrIndexPlaceHolder_iterator operator-(int i)
  {
    VarConstrIndexPlaceHolder_iterator tmp(*this);
    while (i > 0)
    {
      --tmp;
      i--;
    }
    return tmp;
  }

  friend bool operator==(const VarConstrIndexPlaceHolder_iterator& it1,
          const VarConstrIndexPlaceHolder_iterator& it2)
  {
    return it1._vcPhPtr == it2._vcPhPtr;
  }

  friend bool operator!=(const VarConstrIndexPlaceHolder_iterator& it1,
          const VarConstrIndexPlaceHolder_iterator& it2)
  {
    return it1._vcPhPtr != it2._vcPhPtr;
  }

} ;

/**
 * Describe a sublist of VarConstrIndexPlaceHolder.
 */
template<typename VCTYPE>
class VarConstrIndexSublist
{
private:
  VarConstrIndexPlaceHolder<VCTYPE> * _head;
  VarConstrIndexPlaceHolder<VCTYPE> * _tail;
  VcIndexStatus::VcStatus _vcIndexStatus;
  size_t _size;

public:
  typedef VarConstrIndexPlaceHolder_iterator<VCTYPE, false> iterator;
  typedef VarConstrIndexPlaceHolder_iterator<VCTYPE, true> const_iterator;

  VarConstrIndexSublist(const VcIndexStatus::VcStatus& status) :
  _head(new VarConstrIndexPlaceHolder<VCTYPE>(NULL, -1)),
  _tail(new VarConstrIndexPlaceHolder<VCTYPE>(NULL, -1)),
  _vcIndexStatus(status),
  _size(0)
  {
    clear();
  }

  virtual ~VarConstrIndexSublist()
  {
    delete _head;
    delete _tail;
  }

  /**
   * Is the sublist is empty?
   * @return true if empty, false elsewhere.
   */
  bool empty() const
  {
    return (_head->nextVcInProb() == _tail);
  }

  /**
   * Add a VarConstrIndexPlaceHolder at the end of the sublist and update the tail (and head) of the sublist.
   * @param vcPhPtr
   */
  void addVarConstrPlaceHolder(VarConstrIndexPlaceHolder<VCTYPE> * vcPhPtr)
  {
    /// Set the next element of the VarConstrIndexPlaceHolder to the tail.
    vcPhPtr->nextVcInProb(_tail);
    /// Set the previous element to the previous element of the tail. </br>
    /// If the sublist is empty, then _tail->prevVcInProb() is _head.
    vcPhPtr->prevVcInProb(_tail->prevVcInProb());
    /// Update the status of vcPtr to the status of the sublist.
    vcPhPtr->vcPtr()->vcIndexStatus(_vcIndexStatus);

    /// Update the next element of the previous element of the tail
    /// with the new VarConstrIndexPlaceHolder. </br>
    /// If the sublist is empty, then _tail->prevVcInProb() is _head,
    /// so we update the next element of the _head.
    _tail->prevVcInProb()->nextVcInProb(vcPhPtr);
    /// Update the previous element of the tail with the new
    /// VarConstrIndexPlaceHolder.
    _tail->prevVcInProb(vcPhPtr);
    ++_size;
  }

  /**
   * Move a VarConstrIndexPlaceHolder inside the current sublist.
   * @param vcPhPtr
   */
  void moveVarConstrPlaceHolder(VarConstrIndexPlaceHolder<VCTYPE> * vcPhPtr)
  {
    /// Update the next element of the previous element of vcPhPtr. </br>
    /// If vcPhPtr is the only element of the sublist then
    /// vcPhPtr->prevVcInProb() is _head and vcPhPtr->nextVcInProb() is _tail.
    vcPhPtr->prevVcInProb()->nextVcInProb(vcPhPtr->nextVcInProb());
    vcPhPtr->nextVcInProb()->prevVcInProb(vcPhPtr->prevVcInProb());


    /// Add vcPhPtr to the sublist.
    addVarConstrPlaceHolder(vcPhPtr);
  }

  /**
   * Decrement the size of the current sublist.
   */
  void decrSize()
  {
    --_size;
  }

  /**
   * @return Size of the sublist
   */
  size_t size() const
  {
    return _size;
  }

  iterator begin()
  {
    return iterator(_head->nextVcInProb());
  }

  const_iterator begin() const
  {
    return const_iterator(_head->nextVcInProb());
  }

  iterator end()
  {
    return iterator(_tail);
  }

  const_iterator end() const
  {
    return const_iterator(_tail);
  }

  /**
   * @return Sublist head.
   */
  VarConstrIndexPlaceHolder<VCTYPE> * head() const
  {
    return _head;
  }

  /**
   * @return Sublist tail.
   */
  VarConstrIndexPlaceHolder<VCTYPE> * tail() const
  {
    return _tail;
  }

  /**
   * Clear the sublist:
   *    * reset the next and previous VarConstrIndexPlaceHolder
   *      for sublist head and tail.
   *    * reset sublist size.
   */
  void clear()
  {
    _head->nextVcInProb(_tail);
    _head->prevVcInProb(_head);

    _tail->prevVcInProb(_head);
    _tail->nextVcInProb(_tail);
    _size = 0;
  }
} ;

/**
 * This object permits to have a constant access time to data from an id.
 */
template<typename VCTYPE>
class VarConstrIndexManager
{
public:

  typedef VarConstrIndexPlaceHolder_iterator<VCTYPE, false> iterator;
  typedef VarConstrIndexPlaceHolder_iterator<VCTYPE, true> const_iterator;

  VarConstrIndexSublist<VCTYPE> _activeStaticSublist;
  VarConstrIndexSublist<VCTYPE> _inactiveStaticSublist;
  VarConstrIndexSublist<VCTYPE> _unsuitableStaticSublist;
  VarConstrIndexSublist<VCTYPE> _incompatibleStaticSublist;

  VarConstrIndexSublist<VCTYPE> _activeDynamicSublist;
  VarConstrIndexSublist<VCTYPE> _inactiveDynamicSublist;
  VarConstrIndexSublist<VCTYPE> _unsuitableDynamicSublist;
  VarConstrIndexSublist<VCTYPE> _incompatibleDynamicSublist;

  VarConstrIndexSublist<VCTYPE> _activeArtificialSublist;
  VarConstrIndexSublist<VCTYPE> _inactiveArtificialSublist;
  VarConstrIndexSublist<VCTYPE> _unsuitableArtificialSublist;
  VarConstrIndexSublist<VCTYPE> _incompatibleArtificialSublist;

  VarConstrIndexSublist<VCTYPE> _undefinedSublist;

protected:

  std::vector<VarConstrIndexPlaceHolder<VCTYPE> *> _vcPtrVector; /// holds both static nd dynamic varconstr
  std::set<VCTYPE *, LexicographicSorting > _dynamicVcPtrSet; /// holds copy of dynamic varconstr to check previous existence
  bool _useDynamicVarPool;

  /**
   * Push back vcPtr inside the vector of VCTYPE.
   * @param vcPtr
   * @param sublist the sublist to add the vcPtr
   *        (_activeStaticSublist, _inactiveStaticSublist...)
   */
  void pushVCInsideSubList(VCTYPE * vcPtr, VarConstrIndexSublist<VCTYPE>& sublist)
  {
    ///If there is some available VarConstrIndexPlaceHolder.
    if (!_undefinedSublist.empty())
    {
      /// Get the index of the last VarConstrIndexPlaceHolder.
      vcPtr->vcIndexInProb(_undefinedSublist.tail()->prevVcInProb()->vcIndexInProb());
      /// Update vcPtr's VarConstrIndexPlaceHolder with the new vcPtr.
      _vcPtrVector[vcPtr->vcIndexInProb()]->vcPtr(vcPtr);
      /// Move the VarConstrIndexPlaceHolder to the sublist.
      sublist.moveVarConstrPlaceHolder(_vcPtrVector[vcPtr->vcIndexInProb()]);
      _undefinedSublist.decrSize();
    } 
    else /// No VarConstrIndexPlaceHolder available, we need to create a new one.
    {
      vcPtr->vcIndexInProb(_vcPtrVector.size());
      /// Push back a new VarConstrIndexPlaceHolder with the index vcIndexInProb.
      _vcPtrVector.push_back(new VarConstrIndexPlaceHolder<VCTYPE>(vcPtr, vcPtr->vcIndexInProb()));
      _vcPtrVector[vcPtr->vcIndexInProb()]->vcPtr(vcPtr);
      /// Add the VarConstrIndexPlaceHolder inside the sublist.
      sublist.addVarConstrPlaceHolder(_vcPtrVector[vcPtr->vcIndexInProb()]);
    }

  }

  /**
   * Get a sublist from the status and the flag (for the moment, only the status is used).
   * @param status Active, Inactive, Unsuitable or Undefined
   * @param flag 's' (static) or 'd' (dynamic) or 'a' (artificial)
   * @return the corresponding sublist.
   */
  VarConstrIndexSublist<VCTYPE>& getSublistFromStatusAndFlag(const VcIndexStatus::VcStatus & status,
                                                             char flag)
  {
    //std::cout << " getSublistFromStatusAndFlag() flag = " << flag << std::endl;

    switch (flag) 
      {
    case 's':
      switch (status) {
      case VcIndexStatus::Undefined:
        return _undefinedSublist;
      case VcIndexStatus::Active:
        return _activeStaticSublist;
      case VcIndexStatus::Inactive:
        return _inactiveStaticSublist;
      case VcIndexStatus::Unsuitable:
        return _unsuitableStaticSublist;
      case VcIndexStatus::Incompatible:
        return _incompatibleStaticSublist;
      default:
        throw GlobalException("VarConstrIndexManager::getTailSubListFromStatusAndFlag: "
                              "the vcIndexStatus is not supported: " + std::to_string(status));
      }

      break;
    case 'd':
      switch (status) {
      case VcIndexStatus::Undefined:
        return _undefinedSublist;
      case VcIndexStatus::Active:
        return _activeDynamicSublist;
      case VcIndexStatus::Inactive:
        return _inactiveDynamicSublist;
      case VcIndexStatus::Unsuitable:
        return _unsuitableDynamicSublist;
      case VcIndexStatus::Incompatible:
        return _incompatibleDynamicSublist;  
      default:
        throw GlobalException("VarConstrIndexManager::getTailSubListFromStatusAndFlag: "
                              "the vcIndexStatus is not supported: " + std::to_string(status));
      }
      break;

    case 'a':
      switch (status) {
      case VcIndexStatus::Undefined:
        return _undefinedSublist;
      case VcIndexStatus::Active:
        return _activeArtificialSublist;
      case VcIndexStatus::Inactive:
        return _inactiveArtificialSublist;
      case VcIndexStatus::Unsuitable:
        return _unsuitableArtificialSublist;
      case VcIndexStatus::Incompatible:
      return _incompatibleArtificialSublist;  
      default:
        throw GlobalException("VarConstrIndexManager::getTailSubListFromStatusAndFlag: "
                               "the vcIndexStatus is not supported: " + std::to_string(status));
      }
      break;

    default:
      throw GlobalException("VarConstrIndexManager::getTailSubListFromStatusAndFlag: "
                            "this flag is not supported: " + std::to_string(flag));
    }
  }

  /**
   * Get a sublist from the status and the flag (for the moment, only the status is used).
   * @param status Active, Inactive, Unsuitable or Undefined
   * @param flag 's' (static) or 'd' (dynamic) or 'a' (active)
   * @return the corresponding sublist.
   */
  const VarConstrIndexSublist<VCTYPE>& getSublistFromStatusAndFlag(const VcIndexStatus::VcStatus& status, char flag) const
  {
    switch (flag) {
    case 's':
      switch (status) {
      case VcIndexStatus::Undefined:
        return _undefinedSublist;
      case VcIndexStatus::Active:
        return _activeStaticSublist;
      case VcIndexStatus::Inactive:
        return _inactiveStaticSublist;
      case VcIndexStatus::Unsuitable:
        return _unsuitableStaticSublist;
      default:
        throw GlobalException("VarConstrIndexManager::getTailSubListFromStatusAndFlag: "
                              "the vcIndexStatus is not supported: " + std::to_string(status));
      }
      break;

    case 'd':
      switch (status) {
      case VcIndexStatus::Undefined:
        return _undefinedSublist;
      case VcIndexStatus::Active:
        return _activeDynamicSublist;
      case VcIndexStatus::Inactive:
        return _inactiveDynamicSublist;
      case VcIndexStatus::Unsuitable:
        return _unsuitableDynamicSublist;
      default:
        throw GlobalException("VarConstrIndexManager::getTailSubListFromStatusAndFlag: "
                              "the vcIndexStatus is not supported: " + std::to_string(status));
      }
      break;

    case 'a':
      switch (status) {
      case VcIndexStatus::Undefined:
        return _undefinedSublist;
      case VcIndexStatus::Active:
        return _activeArtificialSublist;
      case VcIndexStatus::Inactive:
        return _inactiveArtificialSublist;
      case VcIndexStatus::Unsuitable:
        return _unsuitableArtificialSublist;
      default:
        throw GlobalException("VarConstrIndexManager::getTailSubListFromStatusAndFlag: "
                              "the vcIndexStatus is not supported: " + std::to_string(status));
      }
      break;

    default:
      throw GlobalException("VarConstrIndexManager::getTailSubListFromStatusAndFlag: "
                            "this flag is not supported: " + std::to_string(flag));
    }
  }

  /**
   * Copy a varConstrIndexManager's sublist to the current
   * VarConstrIndexManager's sublist. <br>
   *
   * @param varConstrIndexManager the VarConstrIndexManager's sublist to be
   * copied.
   *
   * @param status the sublist's status
   * (Active, Inactive, Undefined or Unsuitable).
   *
   * @param flag the sublist's flag ('s' or 'd' or 'a').
   */
  inline void copyVarConstrIndexInSublist(const VarConstrIndexManager<VCTYPE>& varConstrIndexManager,
                                          const VcIndexStatus::VcStatus& status,
                                          char flag)
  {
    const_iterator cIt;

    cIt = varConstrIndexManager.begin(status, flag);

    for (; cIt != varConstrIndexManager.end(status, flag); ++cIt)
    {
      insert(*cIt, status);
    }
  }

public:

  VarConstrIndexManager() :
  _activeStaticSublist(VcIndexStatus::Active),
  _inactiveStaticSublist(VcIndexStatus::Inactive),
  _unsuitableStaticSublist(VcIndexStatus::Unsuitable),
  _incompatibleStaticSublist(VcIndexStatus::Incompatible),
  _activeDynamicSublist(VcIndexStatus::Active),
  _inactiveDynamicSublist(VcIndexStatus::Inactive),
  _unsuitableDynamicSublist(VcIndexStatus::Unsuitable),
  _incompatibleDynamicSublist(VcIndexStatus::Incompatible),          
  _activeArtificialSublist(VcIndexStatus::Active),
  _inactiveArtificialSublist(VcIndexStatus::Inactive),
  _unsuitableArtificialSublist(VcIndexStatus::Unsuitable),
  _incompatibleArtificialSublist(VcIndexStatus::Incompatible),        
  _undefinedSublist(VcIndexStatus::Undefined),
  _useDynamicVarPool(false)
  {
  }

  VarConstrIndexManager(const VarConstrIndexManager<VCTYPE>& varConstrIndexManager) :
  _activeStaticSublist(VcIndexStatus::Active),
  _inactiveStaticSublist(VcIndexStatus::Inactive),
  _unsuitableStaticSublist(VcIndexStatus::Unsuitable),
  _incompatibleStaticSublist(VcIndexStatus::Incompatible),
  _activeDynamicSublist(VcIndexStatus::Active),
  _inactiveDynamicSublist(VcIndexStatus::Inactive),
  _unsuitableDynamicSublist(VcIndexStatus::Unsuitable),
  _incompatibleDynamicSublist(VcIndexStatus::Incompatible),          
  _activeArtificialSublist(VcIndexStatus::Active),
  _inactiveArtificialSublist(VcIndexStatus::Inactive),
  _unsuitableArtificialSublist(VcIndexStatus::Unsuitable),
  _incompatibleArtificialSublist(VcIndexStatus::Incompatible),        
  _undefinedSublist(VcIndexStatus::Undefined),
  _useDynamicVarPool(false)
  {
    copy(varConstrIndexManager);
  }

  virtual ~VarConstrIndexManager()
  {
    clean();
  }

  VarConstrIndexManager<VCTYPE>& operator=(const VarConstrIndexManager<VCTYPE>& varConstrIndexManager)
  {
    clean();
    copy(varConstrIndexManager);

    return *this;
  }

  /**
   * Verify if the VarConstrId is valid.
   * @param id VarConstrId
   * @return true if id < _vcPtrVector.size() or id <= -1 </br>
   *         else false;
   */
  bool checkIdValidity(const VcIndexInProb& id) const
  /// either active, inactive, or available
  {
    return checkIdValidity(id.vcIndexInProb());
  }

  /**
   * Verify if the VarConstrId is valid.
   * @param id VarConstrId
   * @return true if id < _vcPtrVector.size() or id <= -1 </br>
   *         else false;
   */
  bool checkIdValidity(const long& id) const
  {
    return ((id >= 0) && (id < static_cast<long> (_vcPtrVector.size())));
  }

  /**
   * Insert the object vcPtr to the vcPtr's id.
   * If the object id is already available then we reactivate the object.
   * @param vcPtr pointer.
   * @param subListType can be Active, Inactive, Unsuitable
   * @return true if inserted, false elsewhere.
   * @throw GlobalException if vcPtr is NULL.
   */
  virtual bool insert(VCTYPE* vcPtr, const VcIndexStatus::VcStatus& status)
  {
    int printlevel = 6;

    if (vcPtr != NULL)
      {
	///Test if vcPtr doesn't have an index yet.
	if (!checkIdValidity(*vcPtr))
	  {
	    pushVCInsideSubList(vcPtr, getSublistFromStatusAndFlag(status, vcPtr->flag()));
	    if (vcPtr->flag() == 'd')
	      {
              if (printL(printlevel))
              {
                  std::cout << "_dynamicVcPtrSet before insertion contains: " << std::endl;
                  for(typename std::set<VCTYPE *, LexicographicSorting >::iterator it = _dynamicVcPtrSet.begin();
                      it != _dynamicVcPtrSet.end(); it++)
                  {
                      std::cout << std::hex << *it << std::dec
                                << " " << getDebugInfo(*it)
                                << std::endl;
                  }
              }

              if (vcPtr->isTypeOf(VcId::MastColumnMask) && ((status == VcIndexStatus::Active) || _useDynamicVarPool))
                  _dynamicVcPtrSet.insert(vcPtr);

              if (printL(printlevel))
              {
                  std::cout << "_dynamicVcPtrSet after insertion contains: " << std::endl;
                  for(typename std::set<VCTYPE *, LexicographicSorting >::iterator it = _dynamicVcPtrSet.begin();
                      it != _dynamicVcPtrSet.end(); it++)
                  {
                      std::cout << std::hex << (long)(*it) << std::dec
                                << " " << getDebugInfo(*it)
                                << std::endl;
                  }
              }
	      }
	  } 
	else /// vcPtr has already an id
	  {
	    if (vcPtr->vcIndexStatus() != status)
	      {
            if (!_useDynamicVarPool && vcPtr->isTypeOf(VcId::MastColumnMask))
              {
	            if (vcPtr->vcIndexStatus() == VcIndexStatus::Active)
                    _dynamicVcPtrSet.erase(vcPtr);
	            else if (status == VcIndexStatus::Active)
                    _dynamicVcPtrSet.insert(vcPtr);
              }
		    /// decrement the size of the previous subList that contains vcPtr;
		    getSublistFromStatusAndFlag(vcPtr->vcIndexStatus(), vcPtr->flag()).decrSize();

		    /// Move vcPtr to another subList.
		    getSublistFromStatusAndFlag(status, vcPtr->flag()).moveVarConstrPlaceHolder(_vcPtrVector[vcPtr->vcIndexInProb()]);
	      }
	  }
      }
    else
    {
        throw GlobalException("VarConstrIndexManager::insert : The VarConstr pointer is NULL");
    }

      return true;
  }

  /**
   * Multiple insert.
   * @param begin
   * @param end
   * @return
   */
  virtual bool insert(iterator begin, iterator end, const VcIndexStatus::VcStatus& status)
  {
    bool res = true;
    for (iterator it = begin; it != end; ++it)
    {
      res &= insert(*it, status);
    }
    return res;
  }


  /**
   * This function is just for compatibility reason with std::set erase function.
   * @param vcPtr
   * @return true if the vcPtr is removed, false otherwise.
   */
  virtual bool erase(VCTYPE* vcPtr)
  {
    if (vcPtr != NULL)
    {
        VcIndexStatus::VcStatus status = vcPtr->vcIndexStatus();
      if (checkIdValidity(vcPtr->vcIndexInProb()))
      {
        /// decrement the size of the previous subList that contains vcPtr;
        getSublistFromStatusAndFlag(vcPtr->vcIndexStatus(), vcPtr->flag()).decrSize();

        _undefinedSublist.moveVarConstrPlaceHolder(_vcPtrVector[vcPtr->vcIndexInProb()]);

      if(vcPtr->flag() == 'd')
      {
        if(printL(7))
        {
          std::cout << "_dynamicVcPtrSet.size() before = " << _dynamicVcPtrSet.size() << std::endl;
          std::cout << "vcPtr removed from _dynamicVcPtrSet "
                  << std::hex << (long) vcPtr << std::dec << std::endl;
        }  
        if(printL(7))
        {
          std::cout << "_dynamicVcPtrSet contains: " << std::endl;
          for(typename std::set<VCTYPE *, LexicographicSorting >::iterator it = _dynamicVcPtrSet.begin();
                  it != _dynamicVcPtrSet.end(); it++)
          {
            std::cout << std::hex << (long) (*it) << std::dec
                    << " " << getDebugInfo(*it)
                    << std::endl;
          }
        }
        if (vcPtr->isTypeOf(VcId::MastColumnMask) && ((status == VcIndexStatus::Active) || _useDynamicVarPool))
              _dynamicVcPtrSet.insert(vcPtr);

        if(printL(7))
          std::cout << "_dynamicVcPtrSet.size() after = " << _dynamicVcPtrSet.size() << std::endl;
      }

        _vcPtrVector[vcPtr->vcIndexInProb()]->vcPtr(NULL);
        vcPtr->vcIndexInProb(-1);

        return true;
      }
    }

    return false;
  }

  /**
   * Remove all elements in the active sublists.
   */
  virtual void clear()
  {
    for (iterator it = begin(VcIndexStatus::Active, 's');
         it != end(VcIndexStatus::Active, 's'); ++it)
    {
      erase(*it);
    }

    for (iterator it = begin(VcIndexStatus::Active, 'd');
         it != end(VcIndexStatus::Active, 'd'); ++it)
    {
      erase(*it);
    }

    for (iterator it = begin(VcIndexStatus::Active, 'a');
         it != end(VcIndexStatus::Active, 'a'); ++it)
    {
      erase(*it);
    }
  }

  /**
   * Return the total size of the vector.
   * @return the total size of the vector.
   */
  size_t size() const
  {
    return _vcPtrVector.size();
  }

  size_t size(const VcIndexStatus::VcStatus& status, char flag) const
  {
    return getSublistFromStatusAndFlag(status, flag).size();
  }

  /**
   * Get the size for the tree sublist with the same status
   * @param status
   * @return the sum of the tree sublist size.
   */
  size_t size(const VcIndexStatus::VcStatus& status) const
  {
    return getSublistFromStatusAndFlag(status, 's').size() +
            getSublistFromStatusAndFlag(status, 'd').size() +
            getSublistFromStatusAndFlag(status, 'a').size();
  }

  /**
   * Count number of times vcPtr is present inside the VarConstIndexManager
   * @param vcPtr
   * @return 0 if not present, >=1 elsewhere.
   */
  virtual size_t count(VCTYPE * vcPtr) const
  {
    if (vcPtr != NULL)
    {
      if (checkIdValidity(*vcPtr))
      {
        //if (_vcPtrVector[vcPtr->vcIndexInProb()]->vcPtr() == vcPtr)
        {
          return 1;
        }
        return 0;
      }

      if (vcPtr->flag() == 'd')
      {
        return _dynamicVcPtrSet.count(vcPtr);
      }
      return 0;
    }
    return 0;
  }

  /**
   * Count number of times vcPtr is present inside the VarConstIndexManager
   * @param vcPtr
   * @param status
   * @return 0 if not present, >=1 elsewhere.
   */
  virtual size_t count(VCTYPE * vcPtr,
                       const VcIndexStatus::VcStatus& status) const
  {
    if (vcPtr != NULL)
    {
      if (checkIdValidity(*vcPtr))
      {
        //if (_vcPtrVector[vcPtr->vcIndexInProb()]->vcPtr() == vcPtr)
        {
          return ((vcPtr->vcIndexStatus() == status) ? 1 : 0);
        }
        return 0;
      }

      if (vcPtr->flag() == 'd')
      {
        typename std::set<VCTYPE *, LexicographicSorting >::const_iterator itset =
                _dynamicVcPtrSet.find(vcPtr);
        if (itset != _dynamicVcPtrSet.end())
        {
          return ((*itset)->vcIndexStatus() == status) ? 1 : 0;
        }
      }
    }
    return 0;
  }

  /**
   * By default, return the first element inside the Active Static VarConstr list.
   * @param type Active, Inactive or Unsuitable
   * @param flag 's' for static VarConst or 'd' for dynamic VarConstr or 'a'
   * for artificial VarConstr.
   * @return iterator of the first element of the list.
   */
  iterator begin(const VcIndexStatus::VcStatus& status, char flag)
  {
    return getSublistFromStatusAndFlag(status, flag).begin();
  }

  /**
   * By default, return the first element inside the Active Static VarConstr list.
   * @param type Active, Inactive or Unsuitable
   * @param flag 's' for static VarConst or 'd' for dynamic VarConstr or 'a'
   * for artificial VarConstr.
   * @return const_iterator of the first element of the list.
   */
  const_iterator begin(const VcIndexStatus::VcStatus& status,
                       char flag) const
  {
    return getSublistFromStatusAndFlag(status, flag).begin();
  }

  /**
   * By default, return Active Static VarConstr tail.
   * @param type Active, Inactive or Unsuitable
   * @param flag 's' for static VarConst or 'd' for dynamic VarConstr or 'a'
   * for artificial VarConstr.
   * @return the tail iterator.
   */
  iterator end(const VcIndexStatus::VcStatus& status, char flag)
  {
    return getSublistFromStatusAndFlag(status, flag).end();
  }

  /**
   * By default, return Active Static VarConstr tail.
   * @param type Active, Inactive or Unsuitable
   * @param flag 's' for static VarConst or 'd' for dynamic VarConstr or 'a'
   * for artificial VarConstr.
   * @return the tail const_iterator.
   */
  const_iterator end(const VcIndexStatus::VcStatus& status,
                     char flag) const
  {
    return getSublistFromStatusAndFlag(status, flag).end();
  }

  iterator end()
  {
    throw GlobalException("Not Implemented Yet");
  }

  const_iterator end() const
  {
    throw GlobalException("Not Implemented Yet");
  }

  /**
   * Operator to get the vctype at the position index.
   * @param index
   * @return
   */
  VCTYPE* operator[](long index) const
  {
    return at(index)->vcPtr();
  }

  VCTYPE* operator[](long index)
  {
    return at(index)->vcPtr();
  }

  /**
   * Return the  VCTYPE * at the position index.
   * @param index
   * @return VarConstrIndexPlaceHolder
   * @throw GlobalException if index is not valid.
   */
  VCTYPE * vcPtr(long index) const
  {
    if (checkIdValidity(index))
    {
      return _vcPtrVector[index].vcPtr();
    }

    throw GlobalException("VarConstrPlaceHolder with the index: " + std::to_string(index) + " not found");
  }

  /**
   * Return the VarConstrIndexPlaceHolder at the position index.
   * @param index
   * @return VarConstrIndexPlaceHolder
   * @throw GlobalException if index is not valid.
   */
  VarConstrIndexPlaceHolder<VCTYPE>* at(long index) const
  {
    if (checkIdValidity(index))
    {
      return _vcPtrVector[index];
    }

    throw GlobalException("VarConstrPlaceHolder with the index: " + std::to_string(index) + " not found");
  }

  /**
   * Return the VarConstrIndexPlaceHolder at the position index.
   * @param index
   * @return VarConstrIndexPlaceHolder
   * @throw GlobalException if index is not valid.
   */
  VarConstrIndexPlaceHolder<VCTYPE>* at(long index)
  {
    if (checkIdValidity(index))
    {
      return _vcPtrVector[index];
    }

    throw GlobalException("VarConstrPlaceHolder with the id: " + std::to_string(index) + " not found");
  }


  void printDynamicVcPtrSet()
  {
      std::cout << "printDynamicVcPtrSet : ";
      for (auto it = _dynamicVcPtrSet.begin(); it != _dynamicVcPtrSet.end(); ++it)
      {
          std::cout << " " << (*it)->name();
      }
      std::cout << std::endl;
  }
  /**
   * Return the pointer vcPtr if its find in the current VarConstrIndexManager or NULL elsewhere.
   * @param vcPtr
   * @return vcPtr if its find, NULL elsewhere.
   */
  VCTYPE* findPtr(VCTYPE* vcPtr)
  {
    if (vcPtr != NULL)
    {
      if (checkIdValidity(vcPtr->vcIndexInProb())) /// its index is defined, hence the vcPtr exist already
      {
        return vcPtr;
      } 
      else if (vcPtr->flag() == 'd') /// its index is not defined, but the vcPtr can be the clone of an existing dynamically generated vcPtr
      {
        typename std::set<VCTYPE *, LexicographicSorting>::iterator itset = _dynamicVcPtrSet.find(vcPtr);

        if (itset != _dynamicVcPtrSet.end())
        {
          return *itset;
        }
      }
    }
    return NULL;
  }


  //  /**
  //   * Return the iterator corresponding to the vcPtr
  //   * @param index
  //   * @return VarConstrIndexPlaceHolder
  //   * @throw GlobalException if vcPtr is NULL.
  //   */
  //  const_iterator find(VCTYPE* vcPtr, const VcIndexStatus::VcStatus& status) const
  //  {
  //    if (vcPtr != NULL)
  //      {
  //	if (checkIdValidity(vcPtr->vcIndexInProb()))
  //	  {
  //	    const_iterator it(begin(status, vcPtr->flag()) + vcPtr->vcIndexInProb());
  //	    return it;
  //	  }
  //	else if (vcPtr->flag() == 'd')
  //	  {
  //	    typename std::set<VCTYPE *, LexicographicSorting>::const_iterator itset =
  //	      _dynamicVcPtrSet.find(vcPtr);
  //	    if (itset != _dynamicVcPtrSet.end())
  //	      {
  //		const_iterator it(
  //				  begin(status, vcPtr->flag()) + (*itset)->vcIndexInProb());
  //		return it;
  //	      }
  //	  }
  //	return end(status, vcPtr->flag());
  //      }
  //
  //    throw GlobalException("VarConstrIndexManager::find: vcPtr is NULL");
  //  }

  /**
   * Get the next element after vcPtr
   * @param vcPtr
   * @return next element after vcPtr.
   * @throw GlobalException if vcPtr is not inside VarConstrIndexManager.
   */
  VCTYPE* next(VCTYPE* vcPtr)
  {
    if (vcPtr != NULL)
    {
      if (checkIdValidity(vcPtr->vcIndexInProb()))
      {
        return _vcPtrVector[vcPtr->vcIndexInProb()]->nextVcInProb()->prevVcInProb()->_vcPtr;
      }
      throw GlobalException(
                            "VarConstrIndexManager::next vcPtr id is not valid.");
    }

    return NULL;
  }

  /**
   * Return if a sublist is empty or not.
   * @param status Active, Inactive, Unsuitable or Undefined.
   * @param flag 's' (static) or 'd' (dynamic).
   * @return true if the sublist is empty, false elsewhere.
   */
  bool empty(const VcIndexStatus::VcStatus& status, char flag) const
  {
    return getSublistFromStatusAndFlag(status, flag).empty();
  }

  /**
   * Remove physically all elements in VarConstrIdManager
   */
  void clean()
  {
    while (!_vcPtrVector.empty())
    {
      if(_vcPtrVector.back()!=NULL)
        delete _vcPtrVector.back();
      _vcPtrVector.pop_back();
    }
    _activeStaticSublist.clear();
    _inactiveStaticSublist.clear();
    _unsuitableStaticSublist.clear();
    _activeDynamicSublist.clear();
    _inactiveDynamicSublist.clear();
    _unsuitableDynamicSublist.clear();
    _undefinedSublist.clear();
  }

  /**
   * Make a copy of varConstrIndexManager to the current VarConstrIndexManager.
   *
   * @param varConstrIndexManager VarConstrIndexManager to be copied.
   */
  void copy(const VarConstrIndexManager<VCTYPE>& varConstrIndexManager)
  {
    copyVarConstrIndexInSublist(varConstrIndexManager, VcIndexStatus::Active, 's');
    copyVarConstrIndexInSublist(varConstrIndexManager, VcIndexStatus::Inactive, 's');
    copyVarConstrIndexInSublist(varConstrIndexManager, VcIndexStatus::Unsuitable, 's');

    copyVarConstrIndexInSublist(varConstrIndexManager, VcIndexStatus::Active, 'd');
    copyVarConstrIndexInSublist(varConstrIndexManager, VcIndexStatus::Inactive, 'd');
    copyVarConstrIndexInSublist(varConstrIndexManager, VcIndexStatus::Unsuitable, 'd');

    copyVarConstrIndexInSublist(varConstrIndexManager, VcIndexStatus::Active, 'a');
    copyVarConstrIndexInSublist(varConstrIndexManager, VcIndexStatus::Inactive, 'a');
    copyVarConstrIndexInSublist(varConstrIndexManager, VcIndexStatus::Unsuitable, 'a');

    _useDynamicVarPool = varConstrIndexManager._useDynamicVarPool;
  }

  template <template <typename, typename> class Container,
  typename Value,
  typename Allocator>
  void copyTo(Container<Value, Allocator>& outputContainer,
              VcIndexStatus::VcStatus status, char flag) const
  {
    const_iterator cIt;
    for (cIt = begin(status, flag); cIt != end(status, flag); ++cIt)
    {
      outputContainer.push_back(*cIt);
    }
  }

  template<typename Sort>
  void copyTo(std::set<VCTYPE *, Sort> & outputSet,
              VcIndexStatus::VcStatus status, char flag) const
  {
    const_iterator cIt;
    for (cIt = begin(status, flag); cIt != end(status, flag); ++cIt)
    {
      outputSet.insert(*cIt);
    }
  }

} ;

#endif /* BCVCINDEXINPROBC_HPP_ */
