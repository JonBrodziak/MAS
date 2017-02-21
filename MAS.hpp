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
            info.ShowData();
//            exit(0);
        }

        inline void Run(variable& f) {
            //            std::cout<<"Calls "<<calls++<<std::endl;
            f = 0.0;
            this->catch_biomass_component = 0.0;
            this->survey_biomass_component = 0.0;
            this->fishery_age_comp_component = 0.0;
            this->survey_age_comp_component = 0.0;


            typename mas::Information<REAL_T>::fleet_iterator fit;

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

            for (int i = 0; i < info.survey_models.size(); i++) {
                info.survey_models[i]->Prepare();
            }

            for (fit = info.fleets.begin(); fit != info.fleets.end(); ++fit) {
                (*fit).second->Prepare();
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

            //            for (ait = info.areas.begin(); ait != info.areas.end(); ait++) {
            //                (*ait).second->ComputeProportions();
            //            }

            for (fit = info.fleets.begin(); fit != info.fleets.end(); ++fit) {
                (*fit).second->ComputeProportions();
            }



            for (int i = 0; i < info.survey_models.size(); i++) {
                info.survey_models[i]->ComputeProportions();
            }

            //            for (int y = 0; y < info.nyears; y++) {
            //                for (int s = 0; s < info.nseasons; s++) {
            //                    for (fit = info.fleets.begin(); fit != info.fleets.end(); ++fit) {
            //                        (*fit).second->EvaluateBiomassComponent(y, s);
            //                        (*fit).second->EvaluateAgeCompComponent(y, s);
            //                    }
            //                    for (int i = 0; i < info.survey_models.size(); i++) {
            //                        info.survey_models[i]->EvaluateBiomassComponent(y, s);
            //                        info.survey_models[i]->EvaluateAgeCompComponent(y, s);
            //                    }
            //                }
            //            }

            for (fit = info.fleets.begin(); fit != info.fleets.end(); ++fit) {
                (*fit).second->ComputeNLLComponents();
                this->catch_biomass_component += (*fit).second->catch_biomass_component;
                this->fishery_age_comp_component -= (*fit).second->fishery_age_comp_component;
            }

            for (int i = 0; i < info.survey_models.size(); i++) {
                info.survey_models[i]->ComputeNLLComponents();
                this->survey_age_comp_component += info.survey_models[i]->survey_age_comp_component;
                this->survey_biomass_component -= info.survey_models[i]->survey_biomass_component;
            }
            f = this->survey_age_comp_component + this->survey_biomass_component + this->fishery_age_comp_component + this->catch_biomass_component;
//            std::cout << "f = " << this->survey_age_comp_component << " + "
//                    << this->survey_biomass_component << " + " << catch_biomass_component << " + " <<
//                    fishery_age_comp_component << "= " << f << std::endl;

// exit(0);
        }

        void Forecast() {

        }

        void Report() {
            std::ofstream out("mas_report.txt");
            out << std::fixed;

            typename mas::Information<REAL_T>::fleet_iterator fit;
            typename mas::Information<REAL_T>::recruitment_model_iterator rit;
            
            typename std::unordered_map<int, std::shared_ptr<mas::Population<double> > >::iterator it;
            std::unordered_map<int, std::shared_ptr<mas::Population<double> > >& pops =
                    info.GetPopulations();
            mas::Information<double>::area_iterator ait;

            
            //prepare the recruitment deviations.
            
            for(rit = info.recruitment_models.begin(); rit != info.recruitment_models.end(); ++rit){
                (*rit).second->Prepare();
            }
            
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

            for (int i = 0; i < info.survey_models.size(); i++) {
                info.survey_models[i]->Prepare();
            }

            for (fit = info.fleets.begin(); fit != info.fleets.end(); ++fit) {
                (*fit).second->Prepare();
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
                out << *(*ait).second;
            }

            for (fit = info.fleets.begin(); fit != info.fleets.end(); ++fit) {
                (*fit).second->ComputeProportions();
                out << *(*fit).second;
            }



            for (int i = 0; i < info.survey_models.size(); i++) {
                info.survey_models[i]->ComputeProportions();
                out << *info.survey_models[i];
            }


            std::ofstream out2("populations.txt");
            out2 << std::fixed;
            for (it = pops.begin(); it != pops.end(); ++it) {
                out2<<*(*it).second<<"\n";
            }
            

        }

    private:

    };


}


#endif /* MAS_MAS_HPP */

