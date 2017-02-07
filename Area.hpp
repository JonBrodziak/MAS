/* 
 * File:   Area.hpp
 * 
 * Author: Matthew Supernaw
 * National Oceanic and Atmospheric Administration
 * National Marine Fisheries Service 
 * Sustainable Fisheries Division
 * St. Petersburg, FL, 33701
 * 
 * Created on September 16, 2016, 12:36 PM
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

#ifndef MAS_AREA_HPP
#define MAS_AREA_HPP
#include "Recruitment.hpp"
#include "Mortality.hpp"
#include "Growth.hpp"
#include "Fecundity.hpp"
#include "Selectivity.hpp"
#include "Fleet.hpp"
#include "Survey.hpp"

namespace mas {

    template<typename REAL_T>
    struct Area {
        typedef typename mas::VariableTrait<REAL_T>::variable variable;
        std::string name;
        int id;

        int years;
        int seasons;
        int ages;

        int growth_model_id;
        int recruitment_model_id;
        int mortality_model_id;
        int fecundity_model_id;


        std::vector<variable> N;
        std::vector<variable> N_total;
        std::vector<variable> N_Proportion;
        std::vector<variable> C_Proportion;
        std::vector<variable> C_Biomass_Proportion;
        std::vector<variable> C;
        std::vector<variable> C_Biomass;
        std::vector<variable> C_Biomass_total;
        std::vector<variable> SN; //survey numbers at age
        std::vector<variable> SN_Biomass; //survey numbers at age
        std::vector<variable> SN_Biomass_total; //survey numbers at age
        std::vector<variable> SN_Proportion; //survey numbers at age
        std::vector<variable> SN_Biomass_Proportion; //survey numbers at age


        std::vector<variable> N_males;
        std::vector<variable> N_total_males;
        std::vector<variable> N_Proportion_males;
        std::vector<variable> C_Proportion_males;
        std::vector<variable> C_Biomass_Proportion_males;
        std::vector<variable> C_males;
        std::vector<variable> C_Biomass_males;
        std::vector<variable> C_Biomass_total_males;
        std::vector<variable> SN_males; //survey numbers at age
        std::vector<variable> SN_Biomass_males; //survey numbers at age
        std::vector<variable> SN_Biomass_total_males; //survey numbers at age
        std::vector<variable> SN_Proportion_males; //survey numbers at age
        std::vector<variable> SN_Biomass_Proportion_males; //survey numbers at age

        std::vector<variable> N_females;
        std::vector<variable> N_total_females;
        std::vector<variable> N_Proportion_females;
        std::vector<variable> C_Proportion_females;
        std::vector<variable> C_Biomass_Proportion_females;
        std::vector<variable> C_females;
        std::vector<variable> C_Biomass_females;
        std::vector<variable> C_Biomass_total_females;
        std::vector<variable> SN_females; //survey numbers at age
        std::vector<variable> SN_Biomass_females; //survey numbers at age
        std::vector<variable> SN_Biomass_total_females; //survey numbers at age
        std::vector<variable> SN_Proportion_females; //survey numbers at age
        std::vector<variable> SN_Biomass_Proportion_females; //survey numbers at age


        std::vector<variable> N_diff2;
        std::vector<variable> N_Proportion_diff2;
        std::vector<variable> C_Proportion_diff2;
        std::vector<variable> C_Biomass_Proportion_diff2;
        std::vector<variable> C_diff2;
        std::vector<variable> C_Biomass_diff2;
        std::vector<variable> SN_diff2; //survey numbers at age
        std::vector<variable> SN_Biomass_diff2; //survey numbers at age
        std::vector<variable> SN_Proportion_diff2; //survey numbers at age
        std::vector<variable> SN_Biomass_Proportion_diff2; //survey numbers at age


        std::shared_ptr<GrowthBase<REAL_T> > growth_model;
        std::shared_ptr<RecruitmentBase<REAL_T> > recruitment_model;
        std::shared_ptr<NaturalMortality<REAL_T> > mortality_model;
        std::shared_ptr<FecundityBase<REAL_T> > fecundity_model;
        std::map<int, std::vector< std::shared_ptr<Fleet<REAL_T> > > > seasonal_fleet_operations;

        typedef typename std::map<int, std::vector< std::shared_ptr<Fleet<REAL_T> > > >::iterator seasonal_fleet_operations_iterator;

        std::map<int, std::map<int, std::shared_ptr<mas::Survey<REAL_T> > > > survey_models; //population, season keyed
        typedef typename std::map<int, std::map<int, std::shared_ptr<mas::Survey<REAL_T> > > > ::iterator survey_model_iterator;
        typedef typename std::map<int, std::shared_ptr<mas::Survey<REAL_T> > >::iterator survey_season_iterator;

        std::vector<std::shared_ptr<DataObject<REAL_T> > > survey_biomass_data;
        std::vector<std::shared_ptr<DataObject<REAL_T> > > survey_proportion_data;
        std::vector<std::shared_ptr<DataObject<REAL_T> > > survey_proportion_length_data;
        std::vector<std::shared_ptr<DataObject<REAL_T> > > survey_mean_size_data;

        std::vector<std::shared_ptr<DataObject<REAL_T> > > catch_biomass_data;
        std::vector<std::shared_ptr<DataObject<REAL_T> > > catch_proportion_data;
        std::vector<std::shared_ptr<DataObject<REAL_T> > > catch_proportion_length_data;
        std::vector<std::shared_ptr<DataObject<REAL_T> > > catch_mean_size_data;

        void Initialize(size_t years, size_t seasons, size_t ages) {
            this->years = years;
            this->seasons = seasons;
            this->ages = ages;


            N.resize(years * seasons * ages);
            C.resize(years * seasons * ages);
            C_Biomass.resize(years * seasons * ages);
            N_Proportion.resize(years * seasons * ages);
            C_Proportion.resize(years * seasons * ages);
            C_Biomass_Proportion.resize(years * seasons * ages);
            SN.resize(years * seasons * ages);
            SN_Biomass.resize(years * seasons * ages);
            SN_Proportion.resize(years * seasons * ages);
            SN_Biomass_Proportion.resize(years * seasons * ages);



            N_males.resize(years * seasons * ages);
            C_males.resize(years * seasons * ages);
            C_Biomass_males.resize(years * seasons * ages);
            N_Proportion_males.resize(years * seasons * ages);
            C_Proportion_males.resize(years * seasons * ages);
            C_Biomass_Proportion_males.resize(years * seasons * ages);
            SN_males.resize(years * seasons * ages);
            SN_Biomass_males.resize(years * seasons * ages);
            SN_Proportion_males.resize(years * seasons * ages);
            SN_Biomass_Proportion_males.resize(years * seasons * ages);


            N_females.resize(years * seasons * ages);
            C_females.resize(years * seasons * ages);
            C_Biomass_females.resize(years * seasons * ages);
            N_Proportion_females.resize(years * seasons * ages);
            C_Proportion_females.resize(years * seasons * ages);
            C_Biomass_Proportion_females.resize(years * seasons * ages);
            SN_females.resize(years * seasons * ages);
            SN_Biomass_females.resize(years * seasons * ages);
            SN_Proportion_females.resize(years * seasons * ages);
            SN_Biomass_Proportion_females.resize(years * seasons * ages);

            N_diff2.resize(years * seasons * ages);
            C_diff2.resize(years * seasons * ages);
            C_Biomass_diff2.resize(years * seasons * ages);
            N_Proportion_diff2.resize(years * seasons * ages);
            C_Proportion_diff2.resize(years * seasons * ages);
            C_Biomass_Proportion_diff2.resize(years * seasons * ages);
            SN_diff2.resize(years * seasons * ages);
            SN_Biomass_diff2.resize(years * seasons * ages);
            SN_Proportion_diff2.resize(years * seasons * ages);
            SN_Biomass_Proportion_diff2.resize(years * seasons * ages);
            this->C_Biomass_total.resize(years * seasons);
            this->SN_Biomass_total.resize(years * seasons);
        }

        void Prepare() {
            for (int i = 0; i < this->N.size(); i++) {
                //                mas::VariableTrait<REAL_T>::SetValue(SN[i], static_cast<REAL_T> (0.0));
                //                mas::VariableTrait<REAL_T>::SetValue(SN_Biomass[i], static_cast<REAL_T> (0.0));
                //                mas::VariableTrait<REAL_T>::SetValue(N[i], static_cast<REAL_T> (0.0));
                //                mas::VariableTrait<REAL_T>::SetValue(C[i], static_cast<REAL_T> (0.0));
                //                mas::VariableTrait<REAL_T>::SetValue(C_Biomass[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(N[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(C[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(C_Biomass[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(N_Proportion[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(C_Proportion[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(C_Biomass_Proportion[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(SN[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(SN_Biomass[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(SN_Proportion[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(SN_Biomass_Proportion[i], static_cast<REAL_T> (0.0));
                
                mas::VariableTrait<REAL_T>::SetValue(N_males[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(C_males[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(C_Biomass_males[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(N_Proportion_males[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(C_Proportion_males[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(C_Biomass_Proportion_males[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(SN_males[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(SN_Biomass_males[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(SN_Proportion_males[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(SN_Biomass_Proportion_males[i], static_cast<REAL_T> (0.0));
                
                mas::VariableTrait<REAL_T>::SetValue(N_females[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(C_females[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(C_Biomass_females[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(N_Proportion_females[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(C_Proportion_females[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(C_Biomass_Proportion_females[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(SN_females[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(SN_Biomass_females[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(SN_Proportion_females[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(SN_Biomass_Proportion_females[i], static_cast<REAL_T> (0.0));
            }
        }

        void PushNumbersAndBiomass(const std::vector<variable>& SN_, //survey numbers at age
                const std::vector<variable>& SN_Biomass_, //survey numbers at age
                const std::vector<variable>& N_,
                const std::vector<variable>& C_,
                const std::vector<variable>& C_Biomass_, int population) {

            for (int i = 0; i < this->N.size(); i++) {
                SN[i] += SN_[i];
                SN_Biomass[i] += SN_Biomass_[i];
                N[i] += N_[i];
                C[i] += C_[i];
                C_Biomass[i] += C_Biomass_[i];
            }

        }

        inline void ComputeProportions() {
            for (int y = 0; y < this->years; y++) {
                for (int s = 0; s < this->seasons; s++) {

                    variable total_sn;
                    variable& total_sn_b = this->SN_Biomass_total[y * this->seasons + s];
                    total_sn_b = static_cast<REAL_T> (0.0);
                    variable total_n;
                    variable total_c;
                    variable& total_c_b = this->C_Biomass_total[y * this->seasons + s];
                    total_c_b = static_cast<REAL_T> (0.0);
                    for (int a = 0; a <this->ages; a++) {
                        total_sn += SN[y * this->seasons * this->ages + (s) * this->ages + a];
                        total_sn_b += SN_Biomass[y * this->seasons * this->ages + (s) * this->ages + a];
                        total_n += N[y * this->seasons * this->ages + (s) * this->ages + a];
                        total_c += C[y * this->seasons * this->ages + (s) * this->ages + a];
                        total_c_b += C_Biomass[y * this->seasons * this->ages + (s) * this->ages + a];
                    }


                    for (int a = 0; a <this->ages; a++) {
                        SN_Proportion[y * this->seasons * this->ages + (s) * this->ages + a] = SN[y * this->seasons * this->ages + (s) * this->ages + a] / total_sn;
                        SN_Biomass_Proportion[y * this->seasons * this->ages + (s) * this->ages + a] = SN_Biomass[y * this->seasons * this->ages + (s) * this->ages + a] / total_sn_b;
                        N_Proportion[y * this->seasons * this->ages + (s) * this->ages + a] = N[y * this->seasons * this->ages + (s) * this->ages + a] / total_n;
                        C_Proportion[y * this->seasons * this->ages + (s) * this->ages + a] = C[y * this->seasons * this->ages + (s) * this->ages + a] / total_c;
                        C_Biomass_Proportion[y * this->seasons * this->ages + (s) * this->ages + a] = C_Biomass[y * this->seasons * this->ages + (s) * this->ages + a] / total_c_b;
                    }


                }
            }
        }

        const variable Compute() {
            variable f;
            variable c_f;
            variable c_p_f;
            variable s_f;
            variable s_p_f;
            for (int y = 0; y < this->years; y++) {
                for (int s = 0; s < this->seasons; s++) {
                    REAL_T temp = static_cast<REAL_T> (0.0);

                    for (int i = 0; i < this->catch_biomass_data.size(); i++) {
                        temp += catch_biomass_data[i]->get(y, s);
                    }
                    c_f += ((this->C_Biomass_total[y * seasons + s] - temp)*(this->C_Biomass_total[y * seasons + s] - temp));

                    temp = static_cast<REAL_T> (0.0);
                    for (int i = 0; i < this->survey_biomass_data.size(); i++) {
                        temp += survey_biomass_data[i]->get(y, s);
                    }
                    s_f += ((this->C_Biomass_total[y * seasons + s] - temp)*(this->C_Biomass_total[y * seasons + s] - temp));

                    temp = static_cast<REAL_T> (0.0);
                    for (int a = 0; a <this->ages; a++) {
                        for (int i = 0; i < this->catch_proportion_data.size(); i++) {
                            c_p_f += (this->C_Proportion[y * this->seasons * this->ages + (s) * this->ages + a] - this->catch_proportion_data[i]->get(y, s, a))*
                                    (this->C_Proportion[y * this->seasons * this->ages + (s) * this->ages + a] - this->catch_proportion_data[i]->get(y, s, a));
                        }
                    }

                    for (int a = 0; a <this->ages; a++) {
                        for (int i = 0; i < this->catch_proportion_data.size(); i++) {
                            s_p_f += (this->SN_Proportion[y * this->seasons * this->ages + (s) * this->ages + a] - this->catch_proportion_data[i]->get(y, s, a))*
                                    (this->SN_Proportion[y * this->seasons * this->ages + (s) * this->ages + a] - this->catch_proportion_data[i]->get(y, s, a));
                        }
                    }

                }
            }
            f = static_cast<REAL_T> (.5) * atl::log(c_f) + static_cast<REAL_T> (.5) * atl::log(s_f) +
                    static_cast<REAL_T> (.5) * atl::log(s_p_f) + static_cast<REAL_T> (.5) * atl::log(c_p_f);

            return f;
        }


    };

    template<typename REAL_T>
    std::ostream& operator<<(std::ostream& out, const mas::Area<REAL_T>& area) {

        out << "Area:\n";
        out << "Name: " << area.name << "\n";
        out << "Id: " << area.id << "\n";
        out << "Growth Model: " << area.growth_model->id << "\n";
        out << "Recruitment Model: " << area.recruitment_model->id << "\n";
        out << "Mortality Model: " << area.mortality_model->id << "\n\n";
        out << "Catch Biomass Objects: " << area.catch_biomass_data.size() << "\n";
        out << "Catch Proportion Objects: " << area.catch_proportion_data.size() << "\n";
        out << "Catch Proportion Length Objects: " << area.catch_proportion_length_data.size() << "\n";
        out << "Catch Mean Size Objects: " << area.catch_mean_size_data.size() << "\n";
        out << "Survey Biomass Objects: " << area.survey_biomass_data.size() << "\n";
        out << "Survey Proportion Objects: " << area.survey_proportion_data.size() << "\n";
        out << "Survey Proportion Length Objects: " << area.survey_proportion_length_data.size() << "\n";
        out << "Survey Mean Size Objects: " << area.survey_mean_size_data.size() << "\n";
        out << "Area " << area.id << "\n";
        out << "Number of Fish at Age:\n";

        for (int a = 0; a < area.ages; a++) {
            for (int y = 0; y < area.years; y++) {
                for (int s = 0; s < area.seasons; s++) {
                    out << area.N[y * area.seasons * area.ages + (s) * area.ages + a] << " ";
                }
            }
            out << "\n";
        }
        out << "\n\n";
        out << "Area " << area.id << "\n";
        out << "Proportion of Fish at Age:\n";

        for (int a = 0; a < area.ages; a++) {
            for (int y = 0; y < area.years; y++) {
                for (int s = 0; s < area.seasons; s++) {
                    out << area.N_Proportion[y * area.seasons * area.ages + (s) * area.ages + a] << " ";
                }
            }
            out << "\n";
        }
        out << "\n\n";
        out << "Area " << area.id << "\n";
        out << "Catch at Age:\n";

        for (int a = 0; a < area.ages; a++) {
            for (int y = 0; y < area.years; y++) {
                for (int s = 0; s < area.seasons; s++) {
                    out << area.C[y * area.seasons * area.ages + (s) * area.ages + a] << " ";
                }
            }
            out << "\n";
        }
        out << "\n\n";
        out << "Area " << area.id << "\n";
        out << "Catch Proportion at Age:\n";

        for (int a = 0; a < area.ages; a++) {
            for (int y = 0; y < area.years; y++) {
                for (int s = 0; s < area.seasons; s++) {
                    out << area.C_Proportion[y * area.seasons * area.ages + (s) * area.ages + a] << " ";
                }
            }
            out << "\n";
        }
        out << "\n\n";
        out << "Area " << area.id << "\n";
        out << "Catch Biomass:\n";

        for (int a = 0; a < area.ages; a++) {
            for (int y = 0; y < area.years; y++) {
                for (int s = 0; s < area.seasons; s++) {
                    out << area.C_Biomass[y * area.seasons * area.ages + (s) * area.ages + a] << " ";
                }
            }
            out << "\n";
        }
        out << "\n\n";
        out << "Area " << area.id << "\n";
        out << "Survey Numbers at Age:\n";

        for (int a = 0; a < area.ages; a++) {
            for (int y = 0; y < area.years; y++) {
                for (int s = 0; s < area.seasons; s++) {
                    out << area.SN[y * area.seasons * area.ages + (s) * area.ages + a] << " ";
                }
            }
            out << "\n";
        }
        out << "\n\n";
        out << "Area " << area.id << "\n";
        out << "Survey Proportion at Age:\n";

        for (int a = 0; a < area.ages; a++) {
            for (int y = 0; y < area.years; y++) {
                for (int s = 0; s < area.seasons; s++) {
                    out << area.SN_Proportion[y * area.seasons * area.ages + (s) * area.ages + a] << " ";
                }
            }
            out << "\n";
        }
        out << "\n\n";

        out << "Survey Biomass:\n";
        out << "Area " << area.id << "\n";
        for (int a = 0; a < area.ages; a++) {
            for (int y = 0; y < area.years; y++) {
                for (int s = 0; s < area.seasons; s++) {
                    out << area.SN_Biomass[y * area.seasons * area.ages + (s) * area.ages + a] << " ";
                }
            }
            out << "\n";
        }
        out << "\n\n";

        return out;
    }


}


#endif /* MAS_AREA_HPP */

