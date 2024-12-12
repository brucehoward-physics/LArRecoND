/**
 *  @file   CheckTime.h
 *
 *  @brief  Header file for the algorithm checking the time at a given point in the running of Pandora execution
 *
 *  $Log: $
 */

// BH -- with thanks to Andy Chappel as I'm copying and editing an algorithm originally used for something else

#ifndef CHECK_TIME_ALGORITHM_H
#define CHECK_TIME_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

// see https://en.cppreference.com/w/cpp/chrono/c/time
#include<ctime>

namespace lar_content
{

  /**
   *  @brief CheckTimeAlgorithm class
   */
  class CheckTimeAlgorithm : public pandora::Algorithm
  {
  public:
    /**
     *  @brief  Default constructor
     */
    CheckTimeAlgorithm();

  private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_auxiliaryInfo; ///< String that prints given text, so you can tell it where you put it in the sequence...
  };

} // namespace lar_content

#endif // #ifndef CHECK_TIME_ALGORITHM_H