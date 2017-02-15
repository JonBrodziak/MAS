/* 
 * File:   Fleet.hpp
 * 
 * Author: Matthew Supernaw
 * National Oceanic and Atmospheric Administration
 * National Marine Fisheries Service 
 * Sustainable Fisheries Division
 * St. Petersburg, FL, 33701
 * 
 * Created on September 16, 2016, 12:34 PM
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

#ifndef MAS_FLEET_HPP
#define MAS_FLEET_HPP

#include "Common.hpp"
#include "Selectivity.hpp"


namespace mas {

    template<typename REAL_T>
    struct Fleet : mas::ModelObject<REAL_T> {
        typedef typename VariableTrait<REAL_T>::variable variable;
        variable f;
        std::string name;


        int years;
        int seasons;
        int ages;

        std::shared_ptr<DataObject<REAL_T> > catch_biomass_data;
        std::shared_ptr<DataObject<REAL_T> > catch_proportion_at_age_data;
        std::shared_ptr<DataObject<REAL_T> > catch_proportion_at_length_data;
        std::shared_ptr<DataObject<REAL_T> > catch_mean_size_at_age_data;

        //area, year X season
        //        std::map<int, std::vector<variable> > f;

        std::unordered_map<int, std::unordered_map<int, int> > season_area_selectivity_ids;
        std::unordered_map<int, std::unordered_map<int, std::shared_ptr<mas::SelectivityBase<REAL_T> > > > season_area_selectivity;

        std::unordered_map<int, std::unordered_map<int, int> > area_season_selectivity_ids;
        std::unordered_map<int, std::unordered_map<int, std::shared_ptr<mas::SelectivityBase<REAL_T> > > > area_season_selectivity;

        std::unordered_map<int, std::unordered_map<int, int> > area_season_fishing_mortality_ids;
        std::unordered_map<int, std::unordered_map<int, std::shared_ptr<mas::FishingMortality<REAL_T> > > > area_season_fishing_mortality;

        std::unordered_map<int, std::unordered_map<int, int> > season_area_fishing_mortality_ids;
        std::unordered_map<int, std::unordered_map<int, std::shared_ptr<mas::FishingMortality<REAL_T> > > > season_area_fishing_mortality;

        typedef typename std::unordered_map<int, std::unordered_map<int, int> >::iterator season_area_selectivity_ids_iterator;
        typedef typename std::unordered_map<int, std::shared_ptr<mas::SelectivityBase<REAL_T> > >::iterator area_sectivity_iterator;
        typedef typename std::unordered_map<int, std::unordered_map<int, int> >::iterator season_area_id_iterator;
        typedef typename std::unordered_map<int, int>::iterator area_id_iteraor;
        typedef typename std::unordered_map<int, int>::iterator season_id_iteraor;
        typedef typename std::unordered_map<int, std::unordered_map<int, std::shared_ptr<mas::SelectivityBase<REAL_T> > > >::iterator season_area_selectivity_iterator;
        typedef typename std::unordered_map<int, std::unordered_map<int, std::shared_ptr<mas::FishingMortality<REAL_T> > > >::iterator season_area_fishing_mortality_iterator;


        std::vector<variable> catch_biomass_total;

        std::vector<variable> catch_biomass_total_males;

        std::vector<variable> catch_biomass_total_females;

        std::vector<variable> numbers_at_age;
        std::vector<variable> numbers_total;
        std::vector<variable> proportion_at_age;
        std::vector<variable> catch_proportion_at_age;
        std::vector<variable> catch_biomass_proportion_at_age;
        std::vector<variable> catch_at_age;
        std::vector<variable> catch_biomass_at_age;


        std::vector<variable> numbers_at_age_males;
        std::vector<variable> numbers_total_males;
        std::vector<variable> proportion_at_age_males;
        std::vector<variable> catch_proportion_at_age_males;
        std::vector<variable> catch_biomass_proportion_at_age_males;
        std::vector<variable> catch_at_age_males;
        std::vector<variable> catch_biomass_at_age_males;

        std::vector<variable> numbers_at_age_females;
        std::vector<variable> numbers_total_females;
        std::vector<variable> proportion_at_age_females;
        std::vector<variable> catch_proportion_at_age_females;
        std::vector<variable> catch_biomass_proportion_at_age_females;
        std::vector<variable> catch_at_age_females;
        std::vector<variable> catch_biomass_at_age_females;


        std::vector<variable> N_diff2;
        std::vector<variable> N_Proportion_diff2;
        std::vector<variable> C_Proportion_diff2;
        std::vector<variable> C_Biomass_Proportion_diff2;
        std::vector<variable> C_diff2;
        std::vector<variable> C_Biomass_diff2;

        variable catch_biomass_component;
        variable fishery_age_comp_component;

        void Initialize(size_t years, size_t seasons, size_t ages) {
            this->years = years;
            this->seasons = seasons;
            this->ages = ages;


            numbers_at_age.resize(years * seasons * ages);
            catch_at_age.resize(years * seasons * ages);
            catch_biomass_at_age.resize(years * seasons * ages);
            proportion_at_age.resize(years * seasons * ages);
            catch_proportion_at_age.resize(years * seasons * ages);
            catch_biomass_proportion_at_age.resize(years * seasons * ages);



            numbers_at_age_males.resize(years * seasons * ages);
            catch_at_age_males.resize(years * seasons * ages);
            catch_biomass_at_age_males.resize(years * seasons * ages);
            proportion_at_age_males.resize(years * seasons * ages);
            catch_proportion_at_age_males.resize(years * seasons * ages);
            catch_biomass_proportion_at_age_males.resize(years * seasons * ages);


            numbers_at_age_females.resize(years * seasons * ages);
            catch_at_age_females.resize(years * seasons * ages);
            catch_biomass_at_age_females.resize(years * seasons * ages);
            proportion_at_age_females.resize(years * seasons * ages);
            catch_proportion_at_age_females.resize(years * seasons * ages);
            catch_biomass_proportion_at_age_females.resize(years * seasons * ages);

            N_diff2.resize(years * seasons * ages);
            C_diff2.resize(years * seasons * ages);
            C_Biomass_diff2.resize(years * seasons * ages);
            N_Proportion_diff2.resize(years * seasons * ages);
            C_Proportion_diff2.resize(years * seasons * ages);
            C_Biomass_Proportion_diff2.resize(years * seasons * ages);
            catch_biomass_total.resize(years * seasons);
            catch_biomass_total_males.resize(years * seasons);
            catch_biomass_total_females.resize(years * seasons);
        }

        void Prepare() {
             this->fishery_age_comp_component = 0.0;
            this->catch_biomass_component = 0.0;
            for (int i = 0; i < this->numbers_at_age.size(); i++) {
                //                mas::VariableTrait<REAL_T>::SetValue(SN[i], static_cast<REAL_T> (0.0));
                //                mas::VariableTrait<REAL_T>::SetValue(SN_Biomass[i], static_cast<REAL_T> (0.0));
                //                mas::VariableTrait<REAL_T>::SetValue(N[i], static_cast<REAL_T> (0.0));
                //                mas::VariableTrait<REAL_T>::SetValue(C[i], static_cast<REAL_T> (0.0));
                //                mas::VariableTrait<REAL_T>::SetValue(C_Biomass[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(numbers_at_age[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(catch_at_age[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(catch_biomass_at_age[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(proportion_at_age[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(catch_proportion_at_age[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(catch_biomass_proportion_at_age[i], static_cast<REAL_T> (0.0));

                mas::VariableTrait<REAL_T>::SetValue(numbers_at_age_males[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(catch_at_age_males[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(catch_biomass_at_age_males[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(proportion_at_age_males[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(catch_proportion_at_age_males[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(catch_biomass_proportion_at_age_males[i], static_cast<REAL_T> (0.0));

                mas::VariableTrait<REAL_T>::SetValue(numbers_at_age_females[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(catch_at_age_females[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(catch_biomass_at_age_females[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(proportion_at_age_females[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(catch_proportion_at_age_females[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(catch_biomass_proportion_at_age_females[i], static_cast<REAL_T> (0.0));
            }


            for (int i = 0; i < this->catch_biomass_total.size(); i++) {
                mas::VariableTrait<REAL_T>::SetValue(catch_biomass_total[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(catch_biomass_total_males[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(catch_biomass_total_females[i], static_cast<REAL_T> (0.0));
            }
        }

        inline void ComputeProportions() {
            for (int y = 0; y < this->years; y++) {
                for (int s = 0; s < this->seasons; s++) {

                    variable total_n;
                    variable total_n_males;
                    variable total_n_females;
                    variable total_c;
                    variable total_c_males;
                    variable total_c_females;
                    variable& total_c_b = this->catch_biomass_total[y * this->seasons + s];
                    variable& total_c_b_males = this->catch_biomass_total[y * this->seasons + s];
                    variable& total_c_b_females = this->catch_biomass_total[y * this->seasons + s];
                    total_c_b = static_cast<REAL_T> (0.0);
                    size_t index = 0;
                    for (int a = 0; a <this->ages; a++) {
                        index = y * this->seasons * this->ages + (s) * this->ages + a;

                        total_n += numbers_at_age[index];
                        total_n_males += numbers_at_age_males[index];
                        total_n_females += numbers_at_age_females[index];

                        total_c += catch_at_age[index];
                        total_c_males += catch_at_age_males[index];
                        total_c_females += catch_at_age_females[index];

                        total_c_b += catch_biomass_at_age[index];
                        total_c_b_males += catch_biomass_at_age_males[index];
                        total_c_b_females += catch_biomass_at_age_females[index];
                    }


                    for (int a = 0; a <this->ages; a++) {
                        index = y * this->seasons * this->ages + (s) * this->ages + a;
                        proportion_at_age[index] = numbers_at_age[index] / total_n;
                        proportion_at_age_males[index] = numbers_at_age_males[index] / total_n_males;
                        proportion_at_age_females[index] = numbers_at_age_females[index] / total_n_females;

                        catch_proportion_at_age[index] = catch_at_age[index] / total_c;
                        catch_proportion_at_age_males[index] = catch_at_age_males[index] / total_c_males;
                        catch_proportion_at_age_females[index] = catch_at_age_females[index] / total_c_females;

                        catch_biomass_proportion_at_age[index] = catch_biomass_at_age[index] / total_c_b;
                        catch_biomass_proportion_at_age_males[index] = catch_biomass_at_age_males[index] / total_c_b_males;
                        catch_biomass_proportion_at_age_females[index] = catch_biomass_at_age_females[index] / total_c_b_females;

                    }


                }
            }
        }

        void ComputeNLLComponents() {

            this->fishery_age_comp_component = 0.0;
            this->catch_biomass_component = 0.0;

            for (int y = 0; y < this->years; y++) {
                for (int s = 0; s < this->seasons; s++) {
                    REAL_T temp = static_cast<REAL_T> (0.0);
                    REAL_T temp_m = static_cast<REAL_T> (0.0);
                    REAL_T temp_f = static_cast<REAL_T> (0.0);
                    temp += catch_biomass_data->get(y, s);


                    //                    std::cout << "temp = " << temp << "\n";
                    this->catch_biomass_component += .5 * atl::pow(std::log(temp + .00001) - atl::log(this->catch_biomass_total[y * seasons + s] + .00001), 2.0) / .05;

                    temp = static_cast<REAL_T> (0.0);
                    for (int a = 0; a <this->ages; a++) {

                        size_t index = y * this->seasons * this->ages + (s) * this->ages + a;
                        this->fishery_age_comp_component -= atl::pow(this->catch_proportion_at_age_data->get(y, s, a)*(std::log(this->catch_proportion_at_age_data->get(y, s, a) + .0001) - atl::log(this->catch_proportion_at_age[index] + .0001)), 2.0);

                    }

                }
            }
        }

        inline void EvaluateBiomassComponent(int year, int season) {

            REAL_T temp = catch_biomass_data->get(year, season);


            if (temp != catch_biomass_data->missing_value) {
                this->catch_biomass_component += .5 * atl::pow(std::log(temp + .00001) -
                        atl::log(this->catch_biomass_total[year * seasons + season] + .00001), 2.0) / .05;
            }
        }

        inline void EvaluateAgeCompComponent(int year, int season) {
            for (int a = 0; a <this->ages; a++) {
                REAL_T temp = this->catch_proportion_at_age_data->get(year, season, a);
                if (temp != this->catch_proportion_at_age_data->missing_value) {
                    size_t index = year * this->seasons * this->ages + (season) * this->ages + a;
                    this->fishery_age_comp_component -= atl::pow(temp * (std::log(temp + .0001) -
                            atl::log(this->catch_proportion_at_age[index] + .0001)), 2.0);
                }

            }
        }


    };

    template<typename REAL_T >
    std::ostream& operator<<(std::ostream& out, const mas::Fleet<REAL_T>& fleet) {

        out << "Fleet:\n";
        out << "Name: " << fleet.name << "\n";
        out << "Id: " << fleet.id << "\n";


        out << "Area " << fleet.id << "\n";
        out << "Catch at Age:\n";

        for (int a = 0; a < fleet.ages; a++) {
            for (int y = 0; y < fleet.years; y++) {
                for (int s = 0; s < fleet.seasons; s++) {
                    out << fleet.catch_at_age[y * fleet.seasons * fleet.ages + (s) * fleet.ages + a] << " ";
                }
            }
            out << "\n";
        }
        out << "\n\n";
        out << "Area " << fleet.id << "\n";
        out << "Catch at Age:\nMales\n";

        for (int a = 0; a < fleet.ages; a++) {
            for (int y = 0; y < fleet.years; y++) {
                for (int s = 0; s < fleet.seasons; s++) {
                    out << fleet.catch_at_age_males[y * fleet.seasons * fleet.ages + (s) * fleet.ages + a] << " ";
                }
            }
            out << "\n";
        }
        out << "\n\n";

        out << "\n\n";
        out << "Area " << fleet.id << "\n";
        out << "Catch at Age:\nFemales\n";

        for (int a = 0; a < fleet.ages; a++) {
            for (int y = 0; y < fleet.years; y++) {
                for (int s = 0; s < fleet.seasons; s++) {
                    out << fleet.catch_at_age_males[y * fleet.seasons * fleet.ages + (s) * fleet.ages + a] << " ";
                }
            }
            out << "\n";
        }
        out << "\n\n";

        out << "\n\n";
        out << "Area " << fleet.id << "\n";
        out << "Catch Proportion at Age:\nMales and Females\n";

        for (int a = 0; a < fleet.ages; a++) {
            for (int y = 0; y < fleet.years; y++) {
                for (int s = 0; s < fleet.seasons; s++) {
                    out << fleet.catch_proportion_at_age[y * fleet.seasons * fleet.ages + (s) * fleet.ages + a] << " ";
                }
            }
            out << "\n";
        }
        out << "\n\n";

        out << "\n\n";
        out << "Area " << fleet.id << "\n";
        out << "Catch Proportion at Age:\nMales\n";

        for (int a = 0; a < fleet.ages; a++) {
            for (int y = 0; y < fleet.years; y++) {
                for (int s = 0; s < fleet.seasons; s++) {
                    out << fleet.catch_proportion_at_age_males[y * fleet.seasons * fleet.ages + (s) * fleet.ages + a] << " ";
                }
            }
            out << "\n";
        }
        out << "\n\n";

        out << "\n\n";
        out << "Area " << fleet.id << "\n";
        out << "Catch Proportion at Age:\nFemales\n";

        for (int a = 0; a < fleet.ages; a++) {
            for (int y = 0; y < fleet.years; y++) {
                for (int s = 0; s < fleet.seasons; s++) {
                    out << fleet.catch_proportion_at_age_females[y * fleet.seasons * fleet.ages + (s) * fleet.ages + a] << " ";
                }
            }
            out << "\n";
        }
        out << "\n\n";

        out << "Area " << fleet.id << "\n";
        out << "Catch Biomass:\nMales and Females\n";

        for (int a = 0; a < fleet.ages; a++) {
            for (int y = 0; y < fleet.years; y++) {
                for (int s = 0; s < fleet.seasons; s++) {
                    out << fleet.catch_biomass_at_age[y * fleet.seasons * fleet.ages + (s) * fleet.ages + a] << " ";
                }
            }
            out << "\n";
        }
        out << "\n\n";

        out << "\n\n";
        out << "Area " << fleet.id << "\n";
        out << "Catch Biomass:\nMales\n";

        for (int a = 0; a < fleet.ages; a++) {
            for (int y = 0; y < fleet.years; y++) {
                for (int s = 0; s < fleet.seasons; s++) {
                    out << fleet.catch_biomass_at_age_males[y * fleet.seasons * fleet.ages + (s) * fleet.ages + a] << " ";
                }
            }
            out << "\n";
        }
        out << "\n\n";

        out << "\n\n";
        out << "Area " << fleet.id << "\n";
        out << "Catch Biomass:\nFemales\n";

        for (int a = 0; a < fleet.ages; a++) {
            for (int y = 0; y < fleet.years; y++) {
                for (int s = 0; s < fleet.seasons; s++) {
                    out << fleet.catch_biomass_at_age_females[y * fleet.seasons * fleet.ages + (s) * fleet.ages + a] << " ";
                }
            }
            out << "\n";
        }
        out << "\n\n";



        return out;
    }



}


#endif /* MAS_FLEET_HPP */

