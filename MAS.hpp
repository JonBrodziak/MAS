/* 
 * File:   MAS.hpp
 * 
 * Author: Matthew Supernaw
 * National Oceanic and Atmospheric Administration
 * National Marine Fisheries Service 
 * Sustainable Fisheries Division
 * St. Petersburg, FL, 33701
 * 
 * Created on September 16, 2016, 12:35 PM
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

#ifndef MAS_MAS_HPP
#define MAS_MAS_HPP

#include "Information.hpp"

namespace mas {

    template<typename REAL_T>
    class MAS {
        typedef typename VariableTrait<REAL_T>::variable variable;
        std::unordered_map<int, mas::Population<REAL_T> > populations;



        std::string data_file;
        std::string config_file;

        //totals 
        std::vector<variable> N;
        std::vector<variable> N_Proportions;
        std::vector<variable> C;
        std::vector<variable> C_Biomass;
        std::vector<variable> C_Proportions;
        std::vector<variable> SurveyNumbers;
        std::vector<variable> Survey_Biomass;
        std::vector<variable> Survey_Proportions;
        int calls = 0;
        variable catch_biomass_component;
        variable survey_biomass_component;
        variable fishery_age_comp_component;
        variable survey_age_comp_component;
        variable recruitment_deviations_component;
    public:
        mas::Information< REAL_T> info;

        MAS() {
        }

        MAS(const MAS<REAL_T>& orig) {
        }

        virtual ~MAS() {
        }

        void Initialize(const std::string& config_file, const std::string& data_file) {
            info.ParseConfig(config_file);
            info.ParseData(data_file);
            info.CreateModel();
        }

        inline void Run(variable& f) {
            //            std::cout<<"Calls "<<calls++<<std::endl;
            f = 0.0;
            this->catch_biomass_component = 0.0;
            this->survey_biomass_component = 0.0;
            this->fishery_age_comp_component = 0.0;
            this->survey_age_comp_component = 0.0;


            typename std::unordered_map<int, std::shared_ptr<mas::Population<double> > >::iterator it;
            std::unordered_map<int, std::shared_ptr<mas::Population<double> > >& pops =
                    info.GetPopulations();
            mas::Information<double>::area_iterator ait;
            /**
             * Prepare areas for evaluation. Resets runtime information.
             */
            for (ait = info.areas.begin(); ait != info.areas.end(); ait++) {
                (*ait).second->Prepare();
            }

            /**
             * Prepare Populations for evaluation. Resets runtime 
             * information.
             */
            for (it = pops.begin(); it != pops.end(); ++it) {
                (*it).second->Prepare();
            }


            /**
             * Evaluate each population and push final numbers to 
             * Area objects.
             */
            for (it = pops.begin(); it != pops.end(); ++it) {
                (*it).second->Evaluate();
            }

            /**
             * Loop through each area and compute proportions for catch, surveys,
             * and numbers. 
             */

            for (ait = info.areas.begin(); ait != info.areas.end(); ait++) {
                (*ait).second->ComputeProportions();
            }

            
            
            this->info.
            
            
//            /**
//             * Evaluate the likelihood function. 
//             */
//            variable t;
//            for (ait = info.areas.begin(); ait != info.areas.end(); ait++) {
//                (*ait).second->ComputeProportions();
//                t += (*ait).second->Compute();
//                this->catch_biomass_component += (*ait).second->catch_biomass_component;
//                this->survey_biomass_component += (*ait).second->survey_biomass_component;
//                this->fishery_age_comp_component += (*ait).second->fishery_age_comp_component;
//                this->survey_age_comp_component += (*ait).second->survey_age_comp_component;
//            }
//
//            f = .5 * atl::log(this->catch_biomass_component) +
//                    .5 * atl::log(this->survey_biomass_component) +
//                    .5 * atl::log(this->fishery_age_comp_component) +
//                    .5 * atl::log(this->survey_age_comp_component);


            //            if(f.GetValue() != f.GetValue()) {
            //                //                for (ait = info.areas.begin(); ait != info.areas.end(); ait++) {
            //                //                    std::cout<<*(*ait).second<<"\n";
            //                //            }
            //                exit(0);
            //            }


        }

        void Forecast() {

        }

        void Report() {
            mas::Information<double>::area_iterator ait;
            std::ofstream out("mas_report.txt");
            out<<std::fixed;
            for (ait = info.areas.begin(); ait != info.areas.end(); ait++) {
                out << *(*ait).second;
            }
        }

    private:

    };


}


#endif /* MAS_MAS_HPP */

