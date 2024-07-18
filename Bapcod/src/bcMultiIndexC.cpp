/**
 *
 * This file bcMultiIndexC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcMultiIndexC.hpp"

#include "bcIndexC.hpp"
#include "bcVarConstrC.hpp"

#ifndef _MSC_VER
/// Only available on Unix
#include <unistd.h>
#include <execinfo.h>
#include <cxxabi.h>
#endif //_MSC_VER

/// @todo Port this code from C to C++
void stackTrace()
{
#ifndef _MSC_VER
    /**
     * Only available on Unix,
     * stack is automaticaly trace in MVS 2008 at exit failure
     */
    int nb_frames = 100;
    fprintf(stderr, "stack trace:\n");

    /// Storage array for stack trace address data
    void* adress[nb_frames + 1];

    std::cout << "======= Memory map: ========" << std::endl;
    /// Retrieve current stack addresses
    int nb_address = backtrace(adress, sizeof(adress) / sizeof(void*));

    if (nb_address == 0)
    {
        fprintf(stderr, "  <empty, possibly corrupt>\n");
        return;
    }

    /**
     * Resolve addresses into strings containing "filename(function+address)",
     * this array must be free()-ed
     */
    char** symbols = backtrace_symbols(adress, nb_address);

    /// Allocate string which will be filled with the demangled function name
    size_t nb_readable_symbols = 256;
    char* readable_symbols = (char*) malloc(nb_readable_symbols);

    /**
     * Iterate over the returned symbol lines. skip the first, it is the
     * address of this function.
     */
    for (int i = 1; i < nb_address; i++)
    {
        char *init_symbol = 0;
        char *init_offset = 0;
        char *end_offset = 0;

        /**
         * Find parentheses and +address offset surrounding the mangled name:
         * ./module(function+0x15c) [0x8048a6d]
         */
        for (char *c = symbols[i]; *c; ++c)
        {
            if (*c == '(')
                init_symbol = c;
            else if (*c == '+')
                init_offset = c;
            else if (*c == ')' && init_offset)
            {
                end_offset = c;
                break;
            }
        }

        if (init_symbol && init_offset && end_offset && init_symbol < init_offset)
        {
            *init_symbol++ = '\0';
            *init_offset++ = '\0';
            *end_offset = '\0';

            int status;
            char* demangle_symbols = abi::__cxa_demangle(init_symbol,
                                                         readable_symbols, &nb_readable_symbols, &status);

            if (status == 0)
            {
                readable_symbols = demangle_symbols;
                fprintf(stderr, "  %s : %s+%s\n", symbols[i], readable_symbols,
                        init_offset);
            }
            else
            {
                /**
                 * Demangling failed. Stderrput function name as a C function with
                 * no arguments.
                 */
                fprintf(stderr, "  %s : %s()+%s\n", symbols[i], init_symbol,
                        init_offset);
            }
        }
        else
        {
            /// Couldn't parse the line? print the whole line.
            fprintf(stderr, "  %s\n", symbols[i]);
        }
    }

    free(readable_symbols);
    free(symbols);
#endif //_MSC_VER
}

MultiIndexNames::MultiIndexNames(const char & firstName,
                const char & secondName,
                const char & thirdName,
                const char & fourthName,
                const char & fifthName,
                const char & sixthName,
                const char & seventhName,
                const char & eighthName)
{
  letterArray[0] = firstName;
  letterArray[1] = secondName;
  letterArray[2] = thirdName;
  letterArray[3] = fourthName;
  letterArray[4] = fifthName;
  letterArray[5] = sixthName;
  letterArray[6] = seventhName;
  letterArray[7] = eighthName;
}

// added by Boris
void MultiIndexNames::append(const MultiIndexNames &idx)
{
	int limit=maxPosition-1;
	int pos=0;
	while (letterArray[pos]!='_' && pos <=limit)
		pos++;
	if (pos>limit)
	{
		std::cerr << "ERROR : MultiIndexNames::append : MultiIndex is too large" << std::endl;
		stackTrace();
		exit(1);
	}
	int pos2=0;
	while (idx.letterArray[pos2]!='_' && pos <=limit)
	{
		letterArray[pos]=idx.letterArray[pos2];
		pos++;pos2++;
	}
	if (pos>limit)
	{
		std::cerr << "ERROR : MultiIndexNames::append : MultiIndex is too large" << std::endl;
		stackTrace();
		exit(1);
	}
}

char MultiIndexNames::first() const
{
  return letterArray[0];
}

char MultiIndexNames::second() const
{
  return letterArray[1];
}

char MultiIndexNames::third() const
{
  return letterArray[2];
}

char MultiIndexNames::fourth() const
{
  return letterArray[3];
}

char MultiIndexNames::fifth() const
{
  return letterArray[4];
}

char MultiIndexNames::sixth() const
{
  return letterArray[5];
}

char MultiIndexNames::seventh() const
{
  return letterArray[6];
}

char MultiIndexNames::eighth() const
{
  return letterArray[7];
}

MultiIndex::MultiIndex(const MultiIndex & id):
    endPosition(id.endPosition)
{
   for (int position = 0; position < 8; ++position)
       indexArray[position] = id.indexArray[position];
}

MultiIndex::MultiIndex(const MultiIndex & id, int numIndicesToDelete):
  endPosition(id.endPosition)
{
  if (numIndicesToDelete > endPosition)
    numIndicesToDelete = endPosition;
  endPosition -= numIndicesToDelete;
  for (int position = numIndicesToDelete; position < endPosition + numIndicesToDelete; ++position)
    {
      indexArray[position - numIndicesToDelete] = id.indexArray[position];
    }
  if (endPosition < maxPosition)
	  indexArray[endPosition] = -1;
}

MultiIndex::MultiIndex(): endPosition(0)
{
  indexArray[0] = -1;
  indexArray[1] = -1;
  indexArray[2] = -1;
  indexArray[3] = -1;
  indexArray[4] = -1;
  indexArray[5] = -1;
  indexArray[6] = -1;
  indexArray[7] = -1;
}

MultiIndex::MultiIndex(int firstIndex): endPosition(1)
{
  indexArray[0] = firstIndex;
  indexArray[1] = -1;
  indexArray[2] = -1;
  indexArray[3] = -1;
  indexArray[4] = -1;
  indexArray[5] = -1;
  indexArray[6] = -1;
  indexArray[7] = -1;
}

MultiIndex::MultiIndex(int firstIndex,
                       int secondIndex): endPosition(2)
{
  indexArray[0] = firstIndex;
  indexArray[1] = secondIndex;
  indexArray[2] = -1;
  indexArray[3] = -1;
  indexArray[4] = -1;
  indexArray[5] = -1;
  indexArray[6] = -1;
  indexArray[7] = -1;
}

MultiIndex::MultiIndex(int firstIndex,
                       int secondIndex,
                       int thirdIndex): endPosition(3)
{
  indexArray[0] = firstIndex;
  indexArray[1] = secondIndex;
  indexArray[2] = thirdIndex;
  indexArray[3] = -1;
  indexArray[4] = -1;
  indexArray[5] = -1;
  indexArray[6] = -1;
  indexArray[7] = -1;
}
MultiIndex::MultiIndex(int firstIndex,
                       int secondIndex,
                       int thirdIndex,
                       int fourthIndex): endPosition(4)
{
  indexArray[0] = firstIndex;
  indexArray[1] = secondIndex;
  indexArray[2] = thirdIndex;
  indexArray[3] = fourthIndex;
  indexArray[4] = -1;
  indexArray[5] = -1;
  indexArray[6] = -1;
  indexArray[7] = -1;
}
MultiIndex::MultiIndex(int firstIndex,
                       int secondIndex,
                       int thirdIndex,
                       int fourthIndex,
                       int fifthIndex): endPosition(5)
{
  indexArray[0] = firstIndex;
  indexArray[1] = secondIndex;
  indexArray[2] = thirdIndex;
  indexArray[3] = fourthIndex;
  indexArray[4] = fifthIndex;
  indexArray[5] = -1;
  indexArray[6] = -1;
  indexArray[7] = -1;
}

MultiIndex::MultiIndex(int firstIndex,
                       int secondIndex,
                       int thirdIndex,
                       int fourthIndex,
                       int fifthIndex,
                       int sixthIndex): endPosition(6)
{
  indexArray[0] = firstIndex;
  indexArray[1] = secondIndex;
  indexArray[2] = thirdIndex;
  indexArray[3] = fourthIndex;
  indexArray[4] = fifthIndex;
  indexArray[5] = sixthIndex;
  indexArray[6] = -1;
  indexArray[7] = -1;
}

MultiIndex::MultiIndex(int firstIndex,
                       int secondIndex,
                       int thirdIndex,
                       int fourthIndex,
                       int fifthIndex,
                       int sixthIndex,
                       int seventhRef) : endPosition(7)
{
  indexArray[0] = firstIndex;
  indexArray[1] = secondIndex;
  indexArray[2] = thirdIndex;
  indexArray[3] = fourthIndex;
  indexArray[4] = fifthIndex;
  indexArray[5] = sixthIndex;
  indexArray[6] = seventhRef;
  indexArray[7] = -1;
}

MultiIndex::MultiIndex(int firstIndex,
                       int secondIndex,
                       int thirdIndex,
                       int fourthIndex,
                       int fifthIndex,
                       int sixthIndex,
                       int seventhRef,
                       int eighthRef) : endPosition(8)
{
  indexArray[0] = firstIndex;
  indexArray[1] = secondIndex;
  indexArray[2] = thirdIndex;
  indexArray[3] = fourthIndex;
  indexArray[4] = fifthIndex;
  indexArray[5] = sixthIndex;
  indexArray[6] = seventhRef;
  indexArray[7] = eighthRef;
}

int MultiIndex::index(int position) const
{
  if (position < 0) return (-1);
  if (position >= maxPosition) return (-1);
  return indexArray[position];
}

int MultiIndex::first() const
{
  return indexArray[0];
}

int MultiIndex::second() const
{
  return indexArray[1];
}

int MultiIndex::third() const
{
  return indexArray[2];
}

int MultiIndex::fourth() const
{
  return indexArray[3];
}

int MultiIndex::fifth() const
{
  return indexArray[4];
}

int MultiIndex::sixth() const
{
  return indexArray[5];
}

int MultiIndex::seventh() const
{
  return indexArray[6];
}

int MultiIndex::eighth() const
{
  return indexArray[7];
}


MultiIndex & MultiIndex::operator+=(int index)
{
	if (endPosition < maxPosition)
		indexArray[endPosition++] = index;
	return *this;
}

bool MultiIndex::operator<(const MultiIndex & that) const
{
  int end(endPosition);
  if (that.endPosition < end) end = that.endPosition;

  for (int position = 0; position < end; ++position)
    {
      if (indexArray[position] < that.indexArray[position])
        return(true);
      if (indexArray[position] > that.indexArray[position])
        return(false);
    }
  /// equal for all

  return (endPosition < that.endPosition);
}



bool MultiIndex::operator==(const MultiIndex & that) const
{
  if (endPosition != that.endPosition)
    return false;

	int i=0;
	while (i<endPosition && indexArray[i]==that.indexArray[i])
		i++;
	if (i==endPosition)
		return true;
	else
		return false;
}


const std::string & MultiIndex::appendRef2name(std::string & name, const MultiIndexNames & indexNames) const
{
  for (int position = 0; position < endPosition; ++position)
    {
      name = std::string(name) + indexNames.letterArray[position] + indexArray[position];
    }

   return(name);
}

const std::string & MultiIndex::appendRefWithBrackets2name(std::string & name) const
{
    if (endPosition == 0)
    {
        name = std::string(name) + "[]";
        return name;
    }
    name = std::string(name) + "[" + indexArray[0];

    for (int position = 1; position < endPosition; ++position)
    {
        name = std::string(name) + "," + indexArray[position];
    }
    name = std::string(name) + "]";

    return(name);
}

std::ostream& operator<<(std::ostream& os, const MultiIndex & that)
{
  for (int position = 0; position < that.endPosition; ++position)
    {
      os <<  "_" << that.indexArray[position];
    }
  return(os);
}

std::string MultiIndex::indicesAsString() const
{
  std::string ret = "";
  
  for (int position = 0; position < endPosition ; ++position)
    {
      ret = ret + indexArray[position];
      ret = ret + " " ;
    }
  
  return(ret);
}

// added by Boris
bool MultiIndex::moveToNext(const MultiIndex &minIndex,const MultiIndex &maxIndex)
{
	++indexArray[endPosition-1];
	for (int i=endPosition-1;(i>0) && (indexArray[i]>maxIndex.indexArray[i]);--i)
	{
		indexArray[i]=minIndex.indexArray[i];
		++indexArray[i-1];
	}
	if (indexArray[0]>maxIndex.indexArray[0])
		return false;
	else
		return true;
}

void MultiIndex::append(const MultiIndex &that)
{
	if (this->endPosition + that.endPosition > maxPosition)
	{
		std::cerr << "ERROR : MultiIndex::append : Index is too long" << std::endl;
		stackTrace();
		exit(1);
	}
	for (int i2=0;i2<that.endPosition;i2++,endPosition++)
	{
		indexArray[endPosition]=that.indexArray[i2];
	}
}

MultiIndex MultiIndex::operator+(const MultiIndex &that) const
{
	MultiIndex idx=*this;
	idx.append(that);
	return idx;
}