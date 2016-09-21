/* 
 * File:   Growth.hpp
 * 
 * Author: Matthew Supernaw
 * National Oceanic and Atmospheric Administration
 * National Marine Fisheries Service 
 * Sustainable Fisheries Division
 * St. Petersburg, FL, 33701
 * 
 * Created on September 16, 2016, 12:33 PM
 * 
 * This File is part of the NOAA, National Marine Fisheries Service 
 * Metapopulation Assessment System project.
 * 
 * This software is a "United States Government Work" under the terms of the
 * United States Copyright Act.  It was written as part of the author's official
 * duties as a United States Government employee and thus cannot be copyrighted.
 * This software is freely available to the public for use. The National Oceanic 
 * And Atmospheric Administration and the U.S. Government have not placed any 
 * restriction on its use or reproduction.  Although all reasonable efforts have 
 * been taken to ensure the accuracy and reliability of the software and data, 
 * the National Oceanic And Atmospheric Administration and the U.S. Government 
 * do not and cannot warrant the performance or results that may be obtained by 
 * using this  software or data. The National Oceanic And Atmospheric
 * Administration and the U.S. Government disclaim all warranties, express or 
 * implied, including warranties of performance, merchantability or fitness 
 * for any particular purpose.
 *
 * Please cite the author(s) in any work or product based on this material.
 *
 */

#ifndef MAS_GROWTH_HPP
#define MAS_GROWTH_HPP

#include "Common.hpp"

namespace mas {

    template<typename REAL_T>
    struct GrowthBase : mas::ModelObject<REAL_T> {
        typedef typename VariableTrait<REAL_T>::variable variable;
        variable a_min;
        variable a_max;
      

        virtual const variable Evaluate(const variable& age) = 0;



    };

    template<typename REAL_T>
    struct VonBertalanffy : GrowthBase<REAL_T> {
        typedef typename VariableTrait<REAL_T>::variable variable;
        variable k;
        variable l_inf;

        const variable Evaluate(const variable& age) {
            return l_inf * (static_cast<REAL_T> (1.0) - atl::exp(-k * (age - this->a_min)));
        }
    };

    template<typename REAL_T>
    struct SchnuteCaseI : GrowthBase<REAL_T > {
        typedef typename VariableTrait<REAL_T>::variable variable;
        variable alpha;
        variable beta;
        variable lmin;
        variable lmax;

        const variable Evaluate(const variable& age) {
            return atl::pow((lmin + (lmax - lmin))*
                    ((static_cast<REAL_T> (1.0) - atl::exp(-alpha * (age - this->a_min))) /
                    (static_cast<REAL_T> (1.0) - atl::exp(-alpha * (this->a_max - this->a_min)))),
                    static_cast<REAL_T> (1.0) / beta);
        }
    };

    template<typename REAL_T>
    struct SchnuteCaseII : GrowthBase<REAL_T > {
        typedef typename VariableTrait<REAL_T>::variable variable;
        variable alpha;
        variable lmin;
        variable lmax;

        virtual const variable Evaluate(const variable& age) {
            return lmin * atl::exp(atl::log(lmax / lmin)*
                    ((static_cast<REAL_T> (1.0) - atl::exp(-alpha * (age - this->a_min))) /
                    (static_cast<REAL_T> (1.0) - atl::exp(-alpha * (this->a_max - this->a_min)))));
        }
    };

    template<typename REAL_T>
    struct SchnuteCaseIII : GrowthBase<REAL_T> {
        typedef typename VariableTrait<REAL_T>::variable variable;
        variable alpha;
        variable beta;
        variable lmin;
        variable lmax;

        const variable Evaluate(const variable& age) {
            return atl::pow((lmin + (lmax - lmin))*
                    ((static_cast<REAL_T> (1.0) - (age - this->a_min)) /
                    (static_cast<REAL_T> (1.0) - (this->a_max - this->a_min))),
                    static_cast<REAL_T> (1.0) / beta);
        }
    };

    template<typename REAL_T>
    struct SchnuteCaseIV : GrowthBase<REAL_T> {
        typedef typename VariableTrait<REAL_T>::variable variable;
        variable alpha;
        variable beta;
        variable lmin;
        variable lmax;

        const variable Evaluate(const variable& age) {
            return lmin * atl::exp(atl::log(lmax / lmin)*
                    ((static_cast<REAL_T> (1.0) - (age - this->a_min)) /
                    (static_cast<REAL_T> (1.0) - (this->a_max - this->a_min))));
        }
    };


}


#endif /* MAS_GROWTH_HPP */
