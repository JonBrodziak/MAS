/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Survey.hpp
 * Author: matthewsupernaw
 *
 * Created on January 30, 2017, 3:22 PM
 */

#ifndef SURVEY_HPP
#define SURVEY_HPP

#include "Selectivity.hpp"
#include "Area.hpp"

namespace mas {


    template<typename REAL_T>
    struct Area;

    template<typename REAL_T>
    struct Survey : public mas::ModelObject<REAL_T> {
        int id;
        int years;
        int seasons;
        int ages;
        typedef typename VariableTrait<REAL_T>::variable variable;

        std::vector<variable> survey_biomass_total; //survey numbers at age

        std::vector<variable> survey_biomass_total_males; //survey numbers at age

        std::vector<variable> survey_biomass_total_females; //survey numbers at age


        std::vector<variable> survey_numbers_at_age; //survey numbers at age
        std::vector<variable> survey_biomass_at_age; //survey numbers at age
        std::vector<variable> survey_proportion_at_age; //survey numbers at age
        std::vector<variable> survey_biomass_proportion_at_age; //survey numbers at age



        std::vector<variable> survey_numbers_at_age_males; //survey numbers at age
        std::vector<variable> survey_biomass_at_age_males; //survey numbers at age
        std::vector<variable> survey_proportion_at_age_males; //survey numbers at age
        std::vector<variable> survey_biomass_proportion_at_age_males; //survey numbers at age


        std::vector<variable> survey_numbers_at_age_females; //survey numbers at age
        std::vector<variable> survey_biomass_at_age_females; //survey numbers at age
        std::vector<variable> survey_proportion_at_age_females; //survey numbers at age
        std::vector<variable> survey_biomass_proportion_at_age_females; //survey numbers at age



        std::vector<variable> SN_diff2; //survey numbers at age
        std::vector<variable> SN_Biomass_diff2; //survey numbers at age
        std::vector<variable> SN_Proportion_diff2; //survey numbers at age
        std::vector<variable> SN_Biomass_Proportion_diff2; //survey numbers at age


        std::string name;
        int population;
        int selectivity_model_id;
        std::vector<int> area_ids;
        std::shared_ptr<DataObject<REAL_T> > survey_biomass_data;
        std::shared_ptr<DataObject<REAL_T> > survey_proportion_at_age_data;
        std::shared_ptr<DataObject<REAL_T> > survey_proportion_at_length_data;
        std::shared_ptr<DataObject<REAL_T> > survey_mean_size_at_age_data;

        std::shared_ptr<mas::SelectivityBase<REAL_T> > selectivity;
        std::unordered_map<int, std::unordered_map<int, int> > area_season_selectivity_ids;
        std::unordered_map<int, std::unordered_map<int, std::shared_ptr<mas::SelectivityBase<REAL_T> > > > area_season_selectivity;

        std::unordered_map<int, std::unordered_map<int, int> > season_area_selectivity_ids;
        std::unordered_map<int, std::unordered_map<int, std::shared_ptr<mas::SelectivityBase<REAL_T> > > > season_area_selectivity;


        std::unordered_map<int, std::unordered_map<int, std::shared_ptr<mas::Area< REAL_T> > > > seasonal_areas_of_operation;

        typedef typename std::unordered_map<int, std::unordered_map<int, int> >::iterator season_area_selectivity_ids_iterator;
        typedef typename std::unordered_map<int, std::unordered_map<int, int> >::const_iterator season_area_selectivity_ids_const_iterator;
        typedef typename std::unordered_map<int, std::shared_ptr<mas::SelectivityBase<REAL_T> > >::iterator area_sectivity_iterator;
        typedef typename std::unordered_map<int, std::unordered_map<int, int> >::iterator season_area_id_iterator;
        typedef typename std::unordered_map<int, int>::iterator area_id_iteraor;
        typedef typename std::unordered_map<int, int>::iterator season_id_iteraor;
        typedef typename std::unordered_map<int, std::unordered_map<int, std::shared_ptr<mas::SelectivityBase<REAL_T> > > >::iterator season_area_selectivity_iterator;


        //runtime
        variable survey_biomass_component;
        variable survey_age_comp_component;

        void Initialize(size_t years, size_t seasons, size_t ages) {
            this->years = years;
            this->seasons = seasons;
            this->ages = ages;



            survey_numbers_at_age.resize(years * seasons * ages);
            survey_biomass_at_age.resize(years * seasons * ages);
            survey_proportion_at_age.resize(years * seasons * ages);
            survey_biomass_proportion_at_age.resize(years * seasons * ages);



            survey_numbers_at_age_males.resize(years * seasons * ages);
            survey_biomass_at_age_males.resize(years * seasons * ages);
            survey_proportion_at_age_males.resize(years * seasons * ages);
            survey_biomass_proportion_at_age_males.resize(years * seasons * ages);

            survey_numbers_at_age_females.resize(years * seasons * ages);
            survey_biomass_at_age_females.resize(years * seasons * ages);
            survey_proportion_at_age_females.resize(years * seasons * ages);
            survey_biomass_proportion_at_age_females.resize(years * seasons * ages);


            SN_diff2.resize(years * seasons * ages);
            SN_Biomass_diff2.resize(years * seasons * ages);
            SN_Proportion_diff2.resize(years * seasons * ages);
            SN_Biomass_Proportion_diff2.resize(years * seasons * ages);
            survey_biomass_total.resize(years * seasons);
            survey_biomass_total_males.resize(years * seasons);
            survey_biomass_total_females.resize(years * seasons);
        }

        void Prepare() {
            for (int i = 0; i < this->survey_numbers_at_age.size(); i++) {

                mas::VariableTrait<REAL_T>::SetValue(survey_numbers_at_age[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(survey_biomass_at_age[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(survey_proportion_at_age[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(survey_biomass_proportion_at_age[i], static_cast<REAL_T> (0.0));


                mas::VariableTrait<REAL_T>::SetValue(survey_numbers_at_age_males[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(survey_biomass_at_age_males[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(survey_proportion_at_age_males[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(survey_biomass_at_age_males[i], static_cast<REAL_T> (0.0));

                mas::VariableTrait<REAL_T>::SetValue(survey_numbers_at_age_females[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(survey_biomass_at_age_females[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(survey_proportion_at_age_females[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(survey_biomass_proportion_at_age_females[i], static_cast<REAL_T> (0.0));
            }


            for (int i = 0; i < this->survey_biomass_total.size(); i++) {
                mas::VariableTrait<REAL_T>::SetValue(survey_biomass_total[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(survey_biomass_total_males[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(survey_biomass_total_females[i], static_cast<REAL_T> (0.0));
            }
        }

        inline void ComputeProportions() {
            for (int y = 0; y < this->years; y++) {
                for (int s = 0; s < this->seasons; s++) {

                    variable total_sn;
                    variable total_sn_males;
                    variable total_sn_females;
                    variable& total_sn_b = this->survey_biomass_total[y * this->seasons + s];
                    variable& total_sn_b_males = this->survey_biomass_total_males[y * this->seasons + s];
                    variable& total_sn_b_females = this->survey_biomass_total_males[y * this->seasons + s];
                    total_sn_b = static_cast<REAL_T> (0.0);
                    total_sn_b_males = static_cast<REAL_T> (0.0);
                    total_sn_b_females = static_cast<REAL_T> (0.0);

                    size_t index = 0;
                    for (int a = 0; a <this->ages; a++) {
                        index = y * this->seasons * this->ages + (s) * this->ages + a;
                        total_sn += survey_numbers_at_age[index];
                        total_sn_males += survey_numbers_at_age_males[index];
                        total_sn_females += survey_numbers_at_age_females[index];

                        total_sn_b += survey_biomass_at_age[index];
                        total_sn_b_males += survey_biomass_at_age_males[index];
                        total_sn_b_females += survey_biomass_at_age_females[index];

                    }


                    for (int a = 0; a <this->ages; a++) {
                        index = y * this->seasons * this->ages + (s) * this->ages + a;
                        survey_proportion_at_age[index] = survey_numbers_at_age[index] / total_sn;
                        survey_proportion_at_age_males[index] = survey_numbers_at_age_males[index] / total_sn_males;
                        survey_proportion_at_age_females[index] = survey_numbers_at_age_females[index] / total_sn_females;

                        survey_biomass_proportion_at_age[index] = survey_biomass_at_age[index] / total_sn_b;
                        survey_biomass_proportion_at_age_males[index] = survey_biomass_at_age_males[index] / total_sn_b_males;
                        survey_biomass_proportion_at_age_females[index] = survey_biomass_at_age_females[index] / total_sn_b_females;



                    }


                }
            }
        }

        void ComputeNLLComponents() {
            REAL_T srv_CV = .2;
            REAL_T o = .00001;
            this->survey_biomass_component = static_cast<REAL_T> (0.0);
            this->survey_age_comp_component = static_cast<REAL_T> (0.0);

            for (int y = 0; y < this->years; y++) {
                for (int s = 0; s < this->seasons; s++) {
                    REAL_T temp = this->survey_biomass_data->get(y, s);
                    if (temp != this->survey_biomass_data->missing_value) {
                        this->survey_biomass_component +=
                                (.5 * atl::pow((std::log(temp + o) -
                                atl::log(this->survey_biomass_total[y * seasons + s] + o)
                                + (std::pow(srv_CV, 2.0) / 2.0)) / srv_CV, 2.0));
//                        std::cout<<"survey_biomass_component = "<<survey_biomass_component<<"\n";
                    }
//(0.5 * SQUARE(((atl::log(obs_srv_5_biomass(i) + o) - atl::log(est_srv_5_biomass(i) + o) + SQUARE((obs_srv_5_CV(i) / 2.0))) / obs_srv_5_CV(i))));
//        }
                    variable sum;
                    for (int a = 0; a <this->ages; a++) {
                        size_t index = y * this->seasons * this->ages + (s) * this->ages + a;
                        sum -=
                                this->survey_proportion_at_age[index] * (atl::log(this->survey_proportion_at_age[index] + o)
                                - std::log(this->survey_proportion_at_age_data->get(y, s, a) + o));

                    }

                    this->survey_age_comp_component -= sum;
//                     std::cout<<"survey_age_comp_component = "<<survey_age_comp_component<<"\n";

                }


            }
              
        }

        inline void EvaluateBiomassComponent(int year, int season) {

            REAL_T srv_CV = .2;
            REAL_T o = .00001;
            REAL_T temp = this->survey_biomass_data->get(year, season);
            if (temp != this->survey_biomass_data->missing_value) {
                this->survey_biomass_component +=
                        (.5) * atl::pow((std::log(temp + o) -
                        atl::log(this->survey_biomass_total[year * seasons + season] + o)
                        + (std::pow(srv_CV, 2.0) / 2.0)) / srv_CV, 2.0);
            }
        }

        inline void EvaluateAgeCompComponent(int year, int season) {
            REAL_T o = .00001;
            variable sum;
            for (int a = 0; a <this->ages; a++) {
                size_t index = year * this->seasons * this->ages + (season) * this->ages + a;
                REAL_T temp = this->survey_proportion_at_age_data->get(year, season, a);
                if (temp != this->survey_proportion_at_age_data->missing_value) {
                    sum -=
                            this->survey_proportion_at_age[index] * (atl::log(this->survey_proportion_at_age[index] + o)
                            - std::log(temp + o));
                }
            }

            this->survey_age_comp_component += sum;
        }

    };

    template<typename REAL_T>
    std::ostream& operator<<(std::ostream& out, const Survey<REAL_T>& survey) {

        typename Survey<REAL_T>::season_area_selectivity_ids_const_iterator sait;
        typename std::unordered_map<int, int>::const_iterator iit;

        out << std::fixed;
        out << "Survey\n";
        out << "Name: " << survey.name << "\n";
        out << "Areas: { \n";

        for (sait = survey.season_area_selectivity_ids.begin(); sait != survey.season_area_selectivity_ids.end(); ++sait) {
            out << "{season: " << (*sait).first << ", areas[ ";
            for (iit = (*sait).second.begin(); iit != (*sait).second.end(); ++iit) {
                out << (*iit).first << " ";
            }
            out << "]}\n";
        }

        out << "}\n";

        out << "# Selectivity: " << survey.area_season_selectivity_ids.size() << "\n";

        out << "Name: " << survey.name << "\n";
        out << "Id: " << survey.id << "\n";


        out << "Survey " << survey.id << "\n";
        out << "Survey Numbers at Age:\nMales and Females\n";

        for (int a = 0; a < survey.ages; a++) {
            for (int y = 0; y < survey.years; y++) {
                for (int s = 0; s < survey.seasons; s++) {
                    out << survey.survey_numbers_at_age[y * survey.seasons * survey.ages + (s) * survey.ages + a] << " ";
                }
            }
            out << "\n";
        }
        out << "\n\n";

        out << "\n\n";
        out << "Survey " << survey.id << "\n";
        out << "Survey Numbers at Age:\nMales\n";

        for (int a = 0; a < survey.ages; a++) {
            for (int y = 0; y < survey.years; y++) {
                for (int s = 0; s < survey.seasons; s++) {
                    out << survey.survey_numbers_at_age_males[y * survey.seasons * survey.ages + (s) * survey.ages + a] << " ";
                }
            }
            out << "\n";
        }
        out << "\n\n";

        out << "\n\n";
        out << "Survey " << survey.id << "\n";
        out << "Survey Numbers at Age:\nFemales\n";

        for (int a = 0; a < survey.ages; a++) {
            for (int y = 0; y < survey.years; y++) {
                for (int s = 0; s < survey.seasons; s++) {
                    out << survey.survey_numbers_at_age_females[y * survey.seasons * survey.ages + (s) * survey.ages + a] << " ";
                }
            }
            out << "\n";
        }
        out << "\n\n";

        out << "Survey " << survey.id << "\n";
        out << "Survey Proportion at Age:\nMales and Females\n";

        for (int a = 0; a < survey.ages; a++) {
            for (int y = 0; y < survey.years; y++) {
                for (int s = 0; s < survey.seasons; s++) {
                    out << survey.survey_proportion_at_age[y * survey.seasons * survey.ages + (s) * survey.ages + a] << " ";
                }
            }
            out << "\n";
        }
        out << "\n\n";

        out << "\n\n";
        out << "Survey " << survey.id << "\n";
        out << "Survey Proportion at Age:\nMales\n";

        for (int a = 0; a < survey.ages; a++) {
            for (int y = 0; y < survey.years; y++) {
                for (int s = 0; s < survey.seasons; s++) {
                    out << survey.survey_proportion_at_age_males[y * survey.seasons * survey.ages + (s) * survey.ages + a] << " ";
                }
            }
            out << "\n";
        }
        out << "\n\n";

        out << "\n\n";
        out << "Survey " << survey.id << "\n";
        out << "Survey Proportion at Age:\nFemales\n";

        for (int a = 0; a < survey.ages; a++) {
            for (int y = 0; y < survey.years; y++) {
                for (int s = 0; s < survey.seasons; s++) {
                    out << survey.survey_proportion_at_age_females[y * survey.seasons * survey.ages + (s) * survey.ages + a] << " ";
                }
            }
            out << "\n";
        }
        out << "\n\n";

        out << "Survey Biomass Total:\n Females\n";
        out << "Survey " << survey.id << "\n";
        for (int a = 0; a < survey.ages; a++) {
            for (int y = 0; y < survey.years; y++) {
                for (int s = 0; s < survey.seasons; s++) {
                    out << survey.survey_biomass_at_age_females[y * survey.seasons * survey.ages + (s) * survey.ages + a] << " ";
                }
            }
            out << "\n";
        }
        out << "\n\n";

        out << "Survey Biomass Total:\n Males\n";
        out << "Survey " << survey.id << "\n";
        for (int a = 0; a < survey.ages; a++) {
            for (int y = 0; y < survey.years; y++) {
                for (int s = 0; s < survey.seasons; s++) {
                    out << survey.survey_biomass_at_age_males[y * survey.seasons * survey.ages + (s) * survey.ages + a] << " ";
                }
            }
            out << "\n";
        }
        out << "\n\n";

        out << "Survey Biomass Total:\n";
        out << "Male and Female\n";
        out << "Survey " << survey.id << "\n";
        for (int a = 0; a < survey.ages; a++) {
            for (int y = 0; y < survey.years; y++) {
                for (int s = 0; s < survey.seasons; s++) {
                    out << survey.survey_biomass_at_age[y * survey.seasons * survey.ages + (s) * survey.ages + a] << " ";
                }
            }
            out << "\n";
        }
        out << "\n\n";

        return out;
    }

}


#endif /* SURVEY_HPP */

