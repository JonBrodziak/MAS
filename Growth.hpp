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

        variable alpha_f = 0.000025;
        variable alpha_m = 0.000025;
        variable beta_f = 3.0;
        variable beta_m = 3.0;

        virtual const variable Evaluate(const variable& age) = 0;

       inline const variable getWeight(const int& sex, const variable& length) {
            return sex == 0 ? alpha_f * atl::pow(length, beta_f) : alpha_m * atl::pow(length, beta_m);
        }

        virtual const std::string Name() {
            return "GrowthBase";
        }

    };

    template<typename REAL_T>
    struct VonBertalanffy : GrowthBase<REAL_T> {
        typedef typename VariableTrait<REAL_T>::variable variable;
        variable k;
        variable l_inf;

        const variable Evaluate(const variable& age) {
            return l_inf * (static_cast<REAL_T> (1.0) - atl::exp(-k * (age - this->a_min)));
        }

        virtual const std::string Name() {
            return "Von Bertalanffy";
        }

        virtual std::string ToString() {
            std::stringstream ss;
            ss << "Von Bertalanffy Growth:\n";
            ss << "k = " << k << "\n";
            ss << "l_inf = " << l_inf << "\n";
            return ss.str();
        }
    };

    template<typename REAL_T>
    struct VonBertalanffyModified : GrowthBase<REAL_T> {
        typedef typename VariableTrait<REAL_T>::variable variable;
        variable lmin;
        variable lmax;
        variable l_inf;
        variable c;

        const variable Evaluate(const variable& age) {
            return lmin + (lmax - lmin)*((static_cast<REAL_T> (1.0) -
                    (atl::pow(c, age - this->a_min))) / (static_cast<REAL_T> (1.0) - atl::pow(c, this->a_max - this->a_min)));
        }

        virtual const std::string Name() {
            return "Von Bertalanffy modified";
        }

        virtual std::string ToString() {
            std::stringstream ss;
            ss << "Modified Von Bertalanffy Growth:\n";
            ss << "lmin = " << lmin << "\n";
            ss << "lmax = " << lmax << "\n";
            ss << "c = " << c << "\n";
            ss << "l_inf = " << l_inf << "\n";
            return ss.str();
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

        virtual const std::string Name() {
            return "Schnute Case I";
        }

        virtual std::string ToString() {
            std::stringstream ss;
            ss << "Scnute Case I Growth:\n";
            ss << "alpha = " << alpha << "\n";
            ss << "beta = " << beta << "\n";
            ss << "lmin = " << lmin << "\n";
            ss << "lmax = " << lmax << "\n";
            return ss.str();
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

        virtual const std::string Name() {
            return "Schnute Case II";
        }

        virtual std::string ToString() {
            std::stringstream ss;
            ss << "Scnute Case II Growth:\n";
            ss << "alpha = " << alpha << "\n";
            ss << "lmin = " << lmin << "\n";
            ss << "lmax = " << lmax << "\n";
            return ss.str();
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

        virtual const std::string Name() {
            return "Schnute Case III";
        }

        virtual std::string ToString() {
            std::stringstream ss;
            ss << "Scnute Case III Growth:\n";
            ss << "alpha = " << alpha << "\n";
            ss << "beta = " << beta << "\n";
            ss << "lmin = " << lmin << "\n";
            ss << "lmax = " << lmax << "\n";
            return ss.str();
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

        virtual const std::string Name() {
            return "Schnute Case IV";
        }

        virtual std::string ToString() {
            std::stringstream ss;
            ss << "Scnute Case IV Growth:\n";
            ss << "alpha = " << alpha << "\n";
            ss << "beta = " << beta << "\n";
            ss << "lmin = " << lmin << "\n";
            ss << "lmax = " << lmax << "\n";
            return ss.str();
        }
    };


}


#endif /* MAS_GROWTH_HPP */

