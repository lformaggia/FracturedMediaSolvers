/*!
 *	@file EclipseFileTranslator.hpp
 *	@brief Template functions to parse Eclipse script file
 *
 *	@author Luca Turconi <lturconi@gmail.com>
 *  @date 12-04-2012
 *
 */

#ifndef ECLIPSEFILETRANSLATOR_HPP_
#define ECLIPSEFILETRANSLATOR_HPP_

#include<fstream>
#include<iostream>
#include<string>
#include<cstring>
#include<vector>
#include<cstdlib>
#include <iomanip>

namespace EclipseFile
{
	
	/*!
		@fn readSection
    	This function allows to read data contained in a section of an Eclipse data file.
    	@param fileName The Eclipse data file
    	@param sectionName The Section to be read
    	@param data The vector where read values are saved
    	@param specialChar Set TRUE if the Eclipse data file is generated on windows system.
    	@param numberOfReadValues The number of values to be read
    	@return TRUE -> operation ended correctly
				FALSE -> an error occurred
    */
	template<class T>
		bool readSection(const std::string fileName, const std::string sectionName,
						 std::vector<T> & data, const bool specialChar=1, const int numberOfReadValues=-1 )
	{
		if(numberOfReadValues!=-1)
			data.resize(numberOfReadValues);
		
		int i=0;
		
		char strBuffer[256];
		char strRep[20];
		int rep;
		char strReal[20];
		T real;
		char * strVal=0;
		char * star=0;
		std::string curLine, keyword;
		std::fstream filestr;
		unsigned int pos, end;
		bool intoSection=0;
		unsigned int npos=-1;
		// bool first=1;
		
		std::cout << " Reading... " << std::endl;
		std::cout << " fileName : " << fileName << std::endl;
		std::cout << " sectionName : " << sectionName << std::endl;
		
		filestr.open (fileName.c_str(), std::fstream::in);
	
		if (filestr.is_open())
		{
			std::cout << std::endl << " File successfully opened" << std::endl;
		}
		else
		{
			std::cerr << std::endl << " *** Error: file not opened *** " << std::endl << std::endl;
			return  0;
		}
		
		while(!filestr.eof() && (numberOfReadValues==-1 || i<numberOfReadValues) )
		{
			std::getline(filestr,curLine);
			if(specialChar)
				curLine.erase(curLine.end()-1);
			// remove the last character ("^M") from the string
			// otherwise it creates problems
			// ^M = character used by windows at the end of each line in text files
			
			// remove comments
			pos = curLine.find("--");
			curLine = curLine.substr(0,pos);
			
			
			if( intoSection==1 )
			{
				std::strcpy(strBuffer,curLine.c_str()); // copy in the c_str type buffer
				strVal = strtok (strBuffer," ");
					
				while( strVal!=0 && std::strcmp(strVal,"/") && (i<numberOfReadValues || numberOfReadValues==-1) )
				{
					star = std::strstr(strVal,"*");
				
					if(star!=0)
					{
						std::strncpy (strRep,strVal,star-strVal);
						strRep[star-strVal]='\0';
						std::strcpy (strReal,star+1);
						rep = std::atoi(strRep);
						real = std::atof(strReal);
					
						if(numberOfReadValues==-1)
						{
							for(int j=0; j<rep; ++j)
							{
								data.push_back( real );
								++i;
							}
						}
						else
						{
							for(int j=0; j<rep && i<numberOfReadValues; ++j)
							{
								data[i] = real;
								++i;
							}
						}
					}
					else
					{
						if(numberOfReadValues==-1)
						{
							data.push_back( (T) std::atof(strVal) );
						}
						else
						{
							data[i] = (T) std::atof(strVal);
						}
						++i;
					}

					strVal = std::strtok(0," ");
					star=0;
				}
					
				if( strVal!=0 && !std::strcmp(strVal,"/") )
				{
					filestr.close();
					std::cout << " File closed" << std::endl;
				
					if( numberOfReadValues==-1 || i==numberOfReadValues )
					{
						std::cout << " " << i << " data extracted from file" << std::endl << std::endl;
						return 1;
					}
					if(i<numberOfReadValues)
					{
						std::cerr << std::endl << " *** Error: Only " << i << " data extracted from file" << std::endl;
						std::cerr << std::endl << "            before end of section is reached" << std::endl;
						return 0;
					}
				}
			}
			
			if( intoSection==0 )
			{
				pos = curLine.find(sectionName);
				
				if( pos!=npos )
				{
					end = curLine.find_first_of(" ",pos);
					// isolate the Eclipse keyword
					if(end!=npos)
					{
						keyword = curLine.substr(pos,end-pos);
					}
					else{
						keyword = curLine.substr(pos);
					}
					
					if( !keyword.compare(sectionName) )
						intoSection=1;
				}
			}
			
		} // end while
		
		filestr.close();
		
		std::cout << " File closed" << std::endl;
		
		if(i==numberOfReadValues || numberOfReadValues==-1)
			// (OR numberOfReadValues==-1) can be omitted
			// if numberOfReadValues=-1 and it reaches this point, surely there is a problem
			// in this case i!=-1 because i>0
		{
			std::cout << " " << i << " data extracted from file" << std::endl << std::endl;
			return 1;
		}
		
		std::cerr << std::endl << " *** Error: Only " << i << " data extracted from file" << std::endl;
		std::cerr << "            before end of file is reached" << std::endl;
		return 0;
		
	}

	/*!
		@fn writeSection
    	This function allows to write data in the Eclipse data file format
    	@param values The vector with the values to be written
    	@param sectionName The Section to be written
    	@param fileName The Eclipse data file to be written
    	@param mode The stream opening mode flags
    	@return TRUE -> operation ended correctly
				FALSE -> an error occurred
    */
	template<class T>
		bool writeSection(std::vector<T> & values,  std::string sectionName="", std::string fileName="a.out",
			 std::ios_base::openmode mode = std::ios_base::out | std::ios_base::app)
	{
		std::fstream filestr;
		filestr << std::scientific << std::setprecision(10);
		
		if(mode & std::ios_base::in)
		{
			std::cerr << " *** Error: try to open in input mode in a write function *** " << std::endl;
			return 0;
		}
		
		filestr.open (fileName.c_str(), mode);
	
		if (filestr.is_open())
		{
			std::cout << std::endl << " File: " << fileName << ", successfully opened";
		}
		else
		{
			std::cerr << std::endl << " *** Error: file not opened *** " << std::endl << std::endl;
			return  0;
		}
		
		// Print open mode type
		if(mode & std::ios_base::app)
			std::cout << " in append mode" << std::endl;
		
		if(mode & std::ios_base::trunc)
		{
			std::cout << " in truncate mode" << std::endl;
			std::cout << "    -> All previous content has been discarded" << std::endl;
		}
		
		std::cout << std::endl << " Writing new section... " << std::endl;
		
		filestr << std::endl;
		filestr << sectionName << std::endl;
		
		int rowMaxElem = 6;
		int rowElemCount = 1;
		for(typename std::vector<T>::iterator it = values.begin(); it!=values.end(); ++it, ++rowElemCount)
		{
			filestr << *it << " ";
			if(rowElemCount==rowMaxElem)
			{
				filestr << std::endl;
				rowElemCount=0;
			}
		}
		
		if(rowElemCount!=1) filestr << std::endl;
		filestr << "/" << std::endl;
		
		std::cout << std::endl << " New section: " << sectionName << std::endl
			<< "  has been successfully created in file: " << fileName << std::endl;
		
		return 1;
	}
	

} // namespace EclipseFile
#endif /* ECLIPSEFILETRANSLATOR_HPP_ */


