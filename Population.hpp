/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Population.hpp
 * Author: matthewsupernaw
 *
 * Created on September 14, 2016, 1:44 PM
 */

#ifndef POPULATION_HPP
#define POPULATION_HPP



#include <memory>
#include <vector>
#include <iomanip>
#include <unordered_map>
#include "Area.hpp"

//#include "../AutoDiff_Standalone/AutoDiff/AutoDiff.hpp"
#include "Movement.hpp"
#include "Recruitment.hpp"

namespace mas {



    template<class REAL_T>
    class Population;

    /**
     *Runtime calculated information for a population in a specific area.
     */
    template<class REAL_T>
    struct AreaPopulationInfo {
        typedef typename mas::VariableTrait<REAL_T>::variable variable;
        int id;
        int uid = -9999;
        /*These need to be defined elsewhere*/
        variable spawning_season_offset = .25;
        variable catch_fraction_of_year = .5;
        variable survey_fraction_of_year = .75;
        variable catch_biomass_units = 1000.0;
        variable survey_biomass_units = 1000.0;
        variable A = 0.000025;
        variable B = 3.0;
        variable F = .1;
        /*==================================*/
        bool males = true;
        mas::FishSexType sex;

        Population<REAL_T>* natal_population;

        bool natal_homing = false;
        std::shared_ptr<Area<REAL_T> > area;
        std::shared_ptr<Area<REAL_T> > natal_area;
        std::shared_ptr<GrowthBase<REAL_T> > growth_model;
        std::shared_ptr<mas::NaturalMortality<REAL_T> > natural_mortality_model; //area specific   
        std::map<int, std::shared_ptr<mas::RecruitmentBase<REAL_T> > > recruitment_model; //season area specific  
        typedef typename std::map<int, std::shared_ptr<mas::RecruitmentBase<REAL_T> > >::iterator recruitment_model_iterator;
        //
        int years;
        int seasons;
        std::vector<REAL_T> ages;


        std::vector<REAL_T> maturity_vector;
        std::vector<variable> length_at_season_start;
        std::vector<variable> length_at_spawning;
        std::vector<variable> length_at_catch_time;
        std::vector<variable> length_at_survey_time;

        std::vector<variable> weight_at_season_start;
        std::vector<variable> weight_at_spawning;
        std::vector<variable> weight_at_catch_time;
        std::vector<variable> weight_at_survey_time;
        std::vector<variable> equilibrium_to_survival_at_spawning;
        std::vector<variable> spawning_biomass_at_age;
        std::vector<variable> spawning_biomass;

        std::vector<variable> imigrants;
        std::vector<variable> redistributed_recruits;
        std::vector<variable> emigrants;
        std::vector<variable> growth;
        std::vector<variable> recruitment;
        std::vector<variable> abundance;

        std::vector<variable> M;
        std::vector<variable> initial_numbers;
        //        std::vector<variable> F;
        std::vector<variable> Z;
        std::vector<variable> S;
        std::vector<variable> SN; //survey numbers at age
        std::vector<variable> SN_Biomass; //survey numbers at age
        std::vector<variable> N;
        std::vector<variable> C;
        std::vector<variable> C_Biomass;
        std::vector<variable> predicted_N;
        //        std::vector<variable> S;

        void Initialize() {
            length_at_season_start.resize(this->ages.size() + 1);
            length_at_spawning.resize(this->ages.size());
            length_at_catch_time.resize(this->ages.size());
            length_at_survey_time.resize(this->ages.size());

            weight_at_season_start.resize(this->ages.size() + 1);
            weight_at_spawning.resize(this->ages.size());
            weight_at_catch_time.resize(this->ages.size());
            weight_at_survey_time.resize(this->ages.size());
            equilibrium_to_survival_at_spawning.resize(this->ages.size());
            spawning_biomass_at_age.resize(this->ages.size());
            spawning_biomass.resize(years * seasons);

            recruitment.resize(years * seasons);
            redistributed_recruits.resize(years * seasons);
            abundance.resize(years * seasons);
            initial_numbers.resize(years * seasons);
            //            F.resize(years * seasons * ages.size());
            emigrants.resize((years) * seasons * ages.size());
            imigrants.resize((years) * seasons * ages.size());
            growth.resize(years * seasons * ages.size());
            Z.resize(years * seasons * ages.size());
            S.resize(years * seasons * ages.size());
            SN.resize(years * seasons * ages.size());
            SN_Biomass.resize(years * seasons * ages.size());
            N.resize(years * seasons * ages.size());
            C.resize(years * seasons * ages.size());
            C_Biomass.resize(years * seasons * ages.size());
            predicted_N.resize(years * seasons * ages.size());

        }

        inline void Reset() {
            for (int i = 0; i < N.size(); i++) {
                mas::VariableTrait<REAL_T>::SetValue(emigrants[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(imigrants[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(growth[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(Z[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(S[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(SN[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(SN_Biomass[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(N[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(C[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(C_Biomass[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(predicted_N[i], static_cast<REAL_T> (0.0));

            }

            for (int i = 0; i < recruitment.size(); i++) {
                mas::VariableTrait<REAL_T>::SetValue(recruitment[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(redistributed_recruits[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(abundance[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(initial_numbers[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(spawning_biomass[i], static_cast<REAL_T> (0.0));

            }

            for (int i = 0; i < this->ages.size(); i++) {
                mas::VariableTrait<REAL_T>::SetValue(length_at_spawning[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(length_at_catch_time[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(length_at_survey_time[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(weight_at_spawning[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(weight_at_catch_time[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(weight_at_survey_time[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(equilibrium_to_survival_at_spawning[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(spawning_biomass_at_age[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(weight_at_season_start[i], static_cast<REAL_T> (0.0));
                mas::VariableTrait<REAL_T>::SetValue(length_at_season_start[i], static_cast<REAL_T> (0.0));

            }

            mas::VariableTrait<REAL_T>::SetValue(weight_at_season_start[this->ages.size()], static_cast<REAL_T> (0.0));
            mas::VariableTrait<REAL_T>::SetValue(length_at_season_start[this->ages.size()], static_cast<REAL_T> (0.0));

        }

        void InitNumbers() {
            for (int a = 0; a < ages.size(); a++) {
                this->N[a] = this->initial_numbers[a];
            }
        }

        inline void IncrementTime(int& y, int& s) {
            if (s == this->seasons) {
                y += 1;
                s = 1;
            } else {
                s++;
            }
        }

        inline void DecrementTime(int& y, int& s) {
            if (s == 1) {
                y -= 1;
                s = seasons;
            } else {
                s--;
            }
        }

        /**
         * Evaluates spawn and recruitment for all ages in a year and season.
         * 
         * @param year
         * @param season
         */
        inline void Recruitment(int year, int season) {
            int y = year;
            int s = season;
            

            //#warning add compiler hint here
            if (year == 0 && season == 1) {
                //                                for (int a = 0; a< this->ages.size(); a++) {
                this->recruitment[year * seasons + (season - 1)] = this->initial_numbers[0];
                //                                }
            } else {
                //                this->DecrementTime(y, s);
                variable a;
                a = this->spawning_biomass[y * seasons + (s - 1)];
                recruitment_model_iterator rit = this->recruitment_model.find(season);
                if (rit != this->recruitment_model.end()) {
                    this->recruitment[year * seasons + (season - 1)] =
                            (*rit).second->Evaluate(a);
                } else {
                    std::cout << "recruitment model not found!!!\n";
                    exit(0);
                }
                this->N[year * this->seasons * this->ages.size() + (season - 1) * this->ages.size()] = this->recruitment[year * seasons + (season - 1)];
            }

        }

        /**
         * Evaluates growth for all ages in a year and season.
         * 
         * @param year
         * @param season
         */
        inline void Growth(int year, int season) {
            //            growth[year * this->seasons * this->ages.size() + (season - 1) * this->ages.size()] = variable(.01);
            //            for (int a = 1; a< this->ages.size(); a++) {
            //                growth[year * this->seasons * this->ages.size() + (season - 1) * this->ages.size() + a] =
            //                        this->area->growth_model->Evaluate(this->ages[a]);
            //            }
            variable fract = static_cast<REAL_T> (season) / static_cast<REAL_T> (this->seasons);
//            std::cout<<length_at_season_start.size()<<"<<<< ---- "<<std::endl;
            length_at_season_start[0] = .01;//variable(.01);
            weight_at_season_start[0] = this->A * atl::exp(this->B * atl::log(length_at_season_start[0]));
            for (int a = 1; a< this->ages.size(); a++) {
                length_at_season_start[a] =
                        this->area->growth_model->Evaluate(this->ages[a] + fract);

                weight_at_season_start[a] = this->area->growth_model->getWeight(this->males, length_at_season_start[a]); //this->A * atl::exp(this->B * atl::log(length_at_season_start[a]));

                length_at_spawning[a - 1] = (static_cast<REAL_T> (1.0) - this->spawning_season_offset) *
                        length_at_season_start[a - 1] + this->spawning_season_offset * length_at_season_start[a];

                length_at_catch_time[a - 1] = (static_cast<REAL_T> (1.0) - this->catch_fraction_of_year) *
                        length_at_season_start[a - 1] + this->catch_fraction_of_year * length_at_season_start[a];

                length_at_survey_time[a - 1] = (static_cast<REAL_T> (1.0) - this->survey_fraction_of_year) *
                        length_at_season_start[a - 1] + this->survey_fraction_of_year * length_at_season_start[a];

                weight_at_spawning[a - 1] = this->A * atl::exp(this->B * atl::log(length_at_spawning[a - 1]));
                weight_at_catch_time[a - 1] = this->A * atl::exp(this->B * atl::log(length_at_catch_time[a - 1]));
                weight_at_survey_time[a - 1] = this->A * atl::exp(this->B * atl::log(length_at_survey_time[a - 1]));
            }

            int index = ages.size() - 1;
            length_at_season_start[ages.size()] =
                    this->area->growth_model->Evaluate(this->ages[ages.size() - 1] + static_cast<REAL_T> (1.0));
            weight_at_season_start[ages.size()] = this->A * atl::exp(this->B * atl::log(length_at_season_start[ages.size()]));

            length_at_spawning[index] = (static_cast<REAL_T> (1.0) - this->spawning_season_offset) *
                    length_at_season_start[index] + this->spawning_season_offset * length_at_season_start[ages.size()];

            length_at_catch_time[index] = (static_cast<REAL_T> (1.0) - this->catch_fraction_of_year) *
                    length_at_season_start[index] + this->catch_fraction_of_year * length_at_season_start[ages.size()];

            length_at_survey_time[index] = (static_cast<REAL_T> (1.0) - this->survey_fraction_of_year) *
                    length_at_season_start[index] + this->survey_fraction_of_year * length_at_season_start[ages.size()];


            weight_at_spawning[index] = this->A * atl::exp(this->B * atl::log(length_at_spawning[index]));
            weight_at_catch_time[index] = this->A * atl::exp(this->B * atl::log(length_at_catch_time[index]));
            weight_at_survey_time[index] = this->A * atl::exp(this->B * atl::log(length_at_survey_time[index]));



        }

        /**
         * Evaluates mortality for all ages in a year and season.
         * @param year
         * @param season
         */
        inline void Mortality(int year, int season) {
            std::vector< std::shared_ptr<Fleet<REAL_T> > >& fleets = this->area->seasonal_fleet_operations[season];

            for (int a = 0; a< this->ages.size(); a++) {


                variable f_a = static_cast<REAL_T> (0.0);
                for (int f = 0; f < fleets.size(); f++) {

                    f_a += fleets[f]->area_season_fishing_mortality[this->area->id][season]->Evaluate(year, (season - 1)) *
                            fleets[f]->season_area_selectivity[season][this->area->id]->Evaluate(ages[a]);
                }

                Z[year * this->seasons * this->ages.size() + (season - 1) * this->ages.size() + a]
                        = this->natural_mortality_model->Evaluate(a) + f_a;
            }
        }

        /**
         * Evaluates fecundity for all ages in a year and season.
         * 
         * @param year
         * @param season
         */
        inline void Fecundity(int year, int season) {
            if (natal_homing) {
                //use natal area parameters
            } else {
                // use  area parameters
            }
        }

        /**
         * Evaluates numbers at age for all ages in a year and season.
         * 
         * @param year
         * @param season
         */
        inline void NumbersAtAge(int year, int season) {



            if (year == 0 && season == 1) {
                //                for (int a = 0; a < ages.size(); a++) {
                //                    this->N[a] = this->initial_numbers[a];
                //                }
            } else {
                int y = year;
                int s = season;
                this->DecrementTime(y, s);

                for (int a = 1; a < ages.size(); a++) {
                    this->N[year * this->seasons * this->ages.size() + (season - 1) * this->ages.size() + a] =
                            this->N[y * this->seasons * this->ages.size() + (s - 1) * this->ages.size() + a - 1] *
                            std::exp(static_cast<REAL_T> (-1.0) * Z[y * this->seasons * this->ages.size() + (s - 1) * this->ages.size() + a]) -
                            this->emigrants[y * this->seasons * this->ages.size() + (s - 1) * this->ages.size() + a - 1] +
                            this->imigrants[y * this->seasons * this->ages.size() + (s - 1) * this->ages.size() + a - 1];
                }
                this->N[year * this->seasons * this->ages.size() + (season - 1) * this->ages.size()] = this->redistributed_recruits[year * seasons + (season - 1)];

                //move fish
                //calc new numbers
            }
            //move fish
            //calc new numbers




        }

        /**
         * Evaluates catch at age for all ages and season
         * @param year
         * @param season
         */
        inline void CatchAtAge(int year, int season) {
            std::vector< std::shared_ptr<Fleet<REAL_T> > >& fleets = this->area->seasonal_fleet_operations[season];
            for (int a = 0; a < this->ages.size(); a++) {
                variable s;
                for (int f = 0; f < fleets.size(); f++) {
                    s += fleets[f]->season_area_selectivity[season][this->area->id]->Evaluate(ages[a]);
                }

                C[year * this->seasons * this->ages.size() + (season - 1) * this->ages.size() + a] =
                        N[year * this->seasons * this->ages.size() + (season - 1) * this->ages.size() + a]*
                        ((this->F * s))*
                        (1.0 - atl::exp(-1.0 * (Z[year * this->seasons * this->ages.size() + (season - 1) * this->ages.size() + a]))) /
                        Z[year * this->seasons * this->ages.size() + (season - 1) * this->ages.size() + a];

                C_Biomass[year * this->seasons * this->ages.size() + (season - 1) * this->ages.size() + a] =
                        C[year * this->seasons * this->ages.size() + (season - 1) * this->ages.size() + a] * this->weight_at_catch_time[a] / this->catch_biomass_units;
            }
        }

        /**
         * Evaluates survey numbers at age for all ages and season
         * @param year
         * @param season
         */
        inline void SurveyNumbersAtAge(int year, int season) {
            typename mas::Area<REAL_T>::survey_model_iterator it = this->area->survey_models.find(this->id);

            if (it != this->area->survey_models.end()) {
                typename mas::Area<REAL_T>::survey_season_iterator sit = (*it).second.find(season);

                if (sit != (*it).second.end()) {
#warning this needs a rewrite ... mapping 
                    for (int a = 0; a < this->ages.size(); a++) {
                        this->SN[year * this->seasons * this->ages.size() + (season - 1) * this->ages.size() + a] =
                                this->N[year * this->seasons * this->ages.size() + (season - 1) * this->ages.size() + a] * atl::exp(-1.0 * survey_fraction_of_year *
                                this->Z[year * this->seasons * this->ages.size() + (season - 1) * this->ages.size() + a])*
                                (*sit).second->season_area_selectivity[season][this->area->id]->Evaluate(this->ages[a]);

                        SN_Biomass[year * this->seasons * this->ages.size() + (season - 1) * this->ages.size() + a] =
                                SN[year * this->seasons * this->ages.size() + (season - 1) * this->ages.size() + a] * this->weight_at_survey_time[a] / this->survey_biomass_units;

                    }
                    //                    std::cout << "has selectivity for survey!!!\n";
                }
            }


        }

        /**
         * Evaluates catch biomass at age for all ages and season
         * @param year
         * @param season
         */
        inline void CatchBiomassAtAge(int year, int season) {
            variable cb = static_cast<REAL_T> (0.0);
            for (int a = 1; a < ages.size(); a++) {


            }
        }

        /**
         * Evaluates catch biomass at age for all ages and season
         * @param year
         * @param season
         */
        inline void SurveyBiomassAtAge(int year, int season) {
            variable sb = static_cast<REAL_T> (0.0);
            for (int a = 1; a < ages.size(); a++) {


            }

        }

        inline void SpawningBiomass(int year, int season) {
            variable sb = static_cast<REAL_T> (0.0);
            for (int a = 1; a < ages.size(); a++) {
                this->equilibrium_to_survival_at_spawning[a] =
                        std::exp(-1.0 * this->spawning_season_offset * Z[year * this->seasons * this->ages.size() + (season - 1) * this->ages.size() + a]);

                this->spawning_biomass_at_age[a] = this->equilibrium_to_survival_at_spawning[a] *
                        this->weight_at_spawning[a] *
                        this->maturity_vector[a];
                sb += this->spawning_biomass_at_age[a] * this->N[year * this->seasons * this->ages.size() + (season - 1) * this->ages.size() + a];

            }
            this->spawning_biomass[year * seasons + (season - 1)] = sb / this->survey_biomass_units;
        }

        inline void PushToArea() {
            this->area->PushNumbersAndBiomass(this->SN,
                    this->SN_Biomass,
                    this->N,
                    this->C,
                    this->C_Biomass,
                    this->id,
                    sex);
        }

    };

    template<typename REAL_T>
    std::ostream& operator<<(std::ostream& out, mas::AreaPopulationInfo<REAL_T>& pi) {
        out << std::fixed;
        out << std::setprecision(2);
        out << "Population " << pi.natal_population->id << "\n";
        out << "Area " << pi.area->id << "\n";
        out << "Total Mortality at Age (Z)\n";
        if (pi.males) {
            out << "Males\n";
        } else {
            out << "Females\n";
        }
        for (int a = 0; a < pi.ages.size(); a++) {
            for (int y = 0; y < pi.years; y++) {
                for (int s = 0; s < pi.seasons; s++) {
                    out << pi.Z[y * pi.seasons * pi.ages.size() + s * pi.ages.size() + a] << " ";
                }
            }
            out << "\n";
        }
        out << "\n\n";

        out << "Population " << pi.natal_population->id << "\n";
        out << "Area " << pi.area->id << "\n";
        out << "Numbers at Age \n";
        if (pi.males) {
            out << "Males\n";
        } else {
            out << "Females\n";
        }
        for (int a = 0; a < pi.ages.size(); a++) {
            for (int y = 0; y < pi.years; y++) {
                for (int s = 0; s < pi.seasons; s++) {
                    out << pi.N[y * pi.seasons * pi.ages.size() + s * pi.ages.size() + a] << " ";
                }
            }
            out << "\n";
        }
        out << "\n\n";

        out << "Population " << pi.natal_population->id << "\n";
        out << "Area " << pi.area->id << "\n";
        out << "Catch Numbers at Age \n";
        if (pi.males) {
            out << "Males\n";
        } else {
            out << "Females\n";
        }
        for (int a = 0; a < pi.ages.size(); a++) {
            for (int y = 0; y < pi.years; y++) {
                for (int s = 0; s < pi.seasons; s++) {
                    out << pi.C[y * pi.seasons * pi.ages.size() + s * pi.ages.size() + a] << " ";
                }
            }
            out << "\n";
        }
        out << "\n\n";

        out << "Population " << pi.natal_population->id << "\n";
        out << "Area " << pi.area->id << "\n";
        out << "Catch Biomass at Age \n";
        if (pi.males) {
            out << "Males\n";
        } else {
            out << "Females\n";
        }
        for (int a = 0; a < pi.ages.size(); a++) {
            for (int y = 0; y < pi.years; y++) {
                for (int s = 0; s < pi.seasons; s++) {
                    out << pi.C_Biomass[y * pi.seasons * pi.ages.size() + s * pi.ages.size() + a] << " ";
                }
            }
            out << "\n";
        }
        out << "\n\n";

        out << "Population " << pi.natal_population->id << "\n";
        out << "Area " << pi.area->id << "\n";
        out << "Survey Numbers at Age \n";
        if (pi.males) {
            out << "Males\n";
        } else {
            out << "Females\n";
        }
        for (int a = 0; a < pi.ages.size(); a++) {
            for (int y = 0; y < pi.years; y++) {
                for (int s = 0; s < pi.seasons; s++) {
                    out << pi.SN[y * pi.seasons * pi.ages.size() + s * pi.ages.size() + a] << " ";
                }
            }
            out << "\n";
        }
        out << "\n\n";

        out << "Population " << pi.natal_population->id << "\n";
        out << "Area " << pi.area->id << "\n";
        out << "Survey Biomass at Age \n";
        if (pi.males) {
            out << "Males\n";
        } else {
            out << "Females\n";
        }
        for (int a = 0; a < pi.ages.size(); a++) {
            for (int y = 0; y < pi.years; y++) {
                for (int s = 0; s < pi.seasons; s++) {
                    out << pi.SN_Biomass[y * pi.seasons * pi.ages.size() + s * pi.ages.size() + a] << " ";
                }
            }
            out << "\n";
        }
        out << "\n\n";

        out << "Population " << pi.natal_population->id << "\n";
        out << "Area " << pi.area->id << "\n";
        out << "Total Immigrants \n";
        if (pi.males) {
            out << "Males\n";
        } else {
            out << "Females\n";
        }
        for (int a = 0; a < pi.ages.size(); a++) {
            for (int y = 0; y < pi.years; y++) {
                for (int s = 0; s < pi.seasons; s++) {
                    out << pi.imigrants[y * pi.seasons * pi.ages.size() + s * pi.ages.size() + a] << " ";
                }
            }
            out << "\n";
        }
        out << "\n\n";

        out << "Population " << pi.natal_population->id << "\n";
        out << "Area " << pi.area->id << "\n";
        out << "Total Emigrants \n";
        if (pi.males) {
            out << "Males\n";
        } else {
            out << "Females\n";
        }
        for (int a = 0; a < pi.ages.size(); a++) {
            for (int y = 0; y < pi.years; y++) {
                for (int s = 0; s < pi.seasons; s++) {
                    out << pi.emigrants[y * pi.seasons * pi.ages.size() + s * pi.ages.size() + a] << " ";
                }
            }
            out << "\n";
        }
        out << "\n\n";

        //        out << "Population " << pi.natal_population->id << "\n";
        //        out << "Area " << pi.area->id << "\n";
        //        out << "Growth Pattern \n";
        //        if (pi.male_chohorts) {
        //            out << "Males\n";
        //        } else {
        //            out << "Females\n";
        //        }
        //        for (int a = 0; a < pi.ages.size(); a++) {
        //            for (int y = 0; y < pi.years; y++) {
        //                for (int s = 0; s < pi.seasons; s++) {
        //                    out << pi.growth[y * pi.seasons * pi.ages.size() + s * pi.ages.size() + a] << " ";
        //                }
        //            }
        //            out << "\n";
        //        }
        //        out << "\n\n";


        out << "Population " << pi.natal_population->id << "\n";
        out << "Area " << pi.area->id << "\n";
        out << "Spawning Biomass \n";
        if (pi.males) {
            out << "Males\n";
        } else {
            out << "Females\n";
        }
        for (int y = 0; y < pi.years; y++) {
            for (int s = 0; s < pi.seasons; s++) {
                out << pi.spawning_biomass[y * pi.seasons + s] << " ";
            }
            out << "\n";
        }


        out << "\n\n";

        out << "Population " << pi.natal_population->id << "\n";
        out << "Area " << pi.area->id << "\n";
        out << "weight_at_spawning \n";
        if (pi.males) {
            out << "Males\n";
        } else {
            out << "Females\n";
        }
        for (int a = 0; a < pi.ages.size(); a++) {

            out << pi.weight_at_spawning[a] << "\n";
        }
        out << "\n\n";
        out << std::fixed;
        out << "Population " << pi.natal_population->id << "\n";
        out << "Area " << pi.area->id << "\n";
        out << "length_at_season_start \n";
        if (pi.males) {
            out << "Males\n";
        } else {
            out << "Females\n";
        }
        //        for (int y = 0; y < pi.years; y++) {
        for (int s = 0; s < pi.ages.size(); s++) {
            out << pi.length_at_season_start[ s] << " ";
            //            }
            out << "\n";

        }
        out << "\n";


        out << "\n\n";
        out << std::fixed;
        out << "Population " << pi.natal_population->id << "\n";
        out << "Area " << pi.area->id << "\n";
        out << "spawning_biomass_at_age \n";
        if (pi.males) {
            out << "Males\n";
        } else {
            out << "Females\n";
        }
        //        for (int y = 0; y < pi.years; y++) {
        for (int s = 0; s < pi.ages.size(); s++) {
            out << pi.spawning_biomass_at_age[ s] << " ";
            //            }
            out << "\n";

        }
        out << "\n";

        out << "\n\n";
        out << std::fixed;
        out << "Population " << pi.natal_population->id << "\n";
        out << "Area " << pi.area->id << "\n";
        out << "equilibrium_to_survival_at_spawning \n";
        if (pi.males) {
            out << "Males\n";
        } else {
            out << "Females\n";
        }
        //        for (int y = 0; y < pi.years; y++) {
        for (int s = 0; s < pi.ages.size(); s++) {
            out << pi.equilibrium_to_survival_at_spawning[ s] << " ";
            //            }
            out << "\n";

        }
        out << "\n";


        out << "\n\n";
        out << std::fixed;
        out << "Population " << pi.natal_population->id << "\n";
        out << "Area " << pi.area->id << "\n";
        out << "Recruitment \n";
        if (pi.males) {
            out << "Males\n";
        } else {
            out << "Females\n";
        }
        //        for (int y = 0; y < pi.years; y++) {
        for (int s = 0; s < pi.recruitment.size(); s++) {
            out << pi.recruitment[ s] << " ";
            //            }
            out << "\n";

        }
        out << "\n";

        out << "\n\n";
        out << std::fixed;
        out << "Redistributed Recruits " << pi.natal_population->id << "\n";
        out << "Area " << pi.area->id << "\n";
        out << "Recruitment \n";
        if (pi.males) {
            out << "Males\n";
        } else {
            out << "Females\n";
        }
        //        for (int y = 0; y < pi.years; y++) {
        for (int s = 0; s < pi.redistributed_recruits.size(); s++) {
            out << pi.redistributed_recruits[ s] << " ";
            //            }
            out << "\n";

        }
        out << "\n";

        return out;
    }

    template<typename REAL_T>
    struct InitialNumbers {
        typedef typename mas::VariableTrait<REAL_T>::variable variable;
        FishSexType type;
        int area_id;
        std::vector<variable> values;
    };

    template<typename REAL_T>
    class Population : public mas::ModelObject<REAL_T> {
    public:
        /*********************************************
         * Area specific natural mortality           *
         *********************************************/
        std::map<int, int> male_natural_mortality_ids; //area, natural mortality model  
        std::map<int, int> female_natural_mortality_ids; //area, natural mortality model  
        typedef typename std::map<int, int>::iterator male_natural_mortality_ids_iterator;
        typedef typename std::map<int, int>::iterator female_natural_mortality_ids_iterator;

        /*********************************************
         * Area specific recruitment                 *
         *********************************************/
        std::map<int, std::map<int, int> > recruitment_ids; //area, recruitment model  
        typedef typename std::map<int, std::map<int, int> >::iterator recruitment_ids_iterator;
        typedef typename std::map<int, int>::iterator recruitment_season_ids_iterator;


        std::vector<InitialNumbers<REAL_T> > initial_numbers;
        typedef typename mas::VariableTrait<REAL_T>::variable variable;
        std::string name;
        int natal_area_id;
        int movement_model_id;
        bool natal_homing = false;
        bool natal_recruitment = false;
        bool move_fish_before_lh = false;
        int years;
        int seasons;
        int areas;
        int ages;
        int growth_id;
        std::vector<variable> SN; //survey numbers at age
        std::vector<variable> SN_Biomass; //survey numbers at age
        std::vector<variable> C;
        std::vector<variable> C_Biomass;

        std::shared_ptr<Area<REAL_T> > natal_area; //birth area
        std::vector<std::shared_ptr<Area<REAL_T> > > areas_list; //all areas

        //Movement Tracking
        typedef typename std::unordered_map<int, AreaPopulationInfo<REAL_T> >::iterator cohort_iterator;
        std::unordered_map<int, AreaPopulationInfo<REAL_T> > males;
        std::unordered_map<int, AreaPopulationInfo<REAL_T> > females;
        std::unordered_map<int, int > movement_models_ids; //season keyed
        typedef typename std::unordered_map<int, int >::iterator movement_model_id_iterator;

        std::shared_ptr<mas::Movement<REAL_T> > movement_model;

        std::unordered_map<int, std::shared_ptr<mas::Movement<REAL_T> > > movement_models; //year keyed
        typedef typename std::unordered_map<int, std::shared_ptr<mas::Movement<REAL_T> > >::iterator movement_model_iterator;
        //Estimable
        std::vector<std::vector<variable> > movement_coefficients;
        std::vector<variable> initial_population_males;
        std::vector<variable> initial_population_females;


        std::unordered_map<int, std::unordered_map<int, std::vector<REAL_T> > > maturity_models; //area / sex
        typedef typename std::unordered_map<int, std::unordered_map<int, std::vector<REAL_T> > >::iterator maturity_models_iterator;
        //    typedef typename std::unordered_map<std::vector<std::vector<variable> > >::iterator movement_coefficient_iterator;

        Population() {
        }

        Population(int years,
                int seasons,
                int areas,
                const std::shared_ptr<Area<REAL_T> >& natal_area,
                const std::vector<std::shared_ptr<Area<REAL_T> > >& areas_list) :
        years(years),
        seasons(seasons),
        areas(areas),
        natal_area(natal_area),
        areas_list(areas) {

            for (int a = 0; a < areas_list.size(); a++) {
                males[areas_list[a]->id].natal_homing = this->natal_homing;
                males[areas_list[a]->id].area = areas_list[a];
                males[areas_list[a]->id].natal_area = this->natal_area;
                males[areas_list[a]->id].Initialize();
                
                females[areas_list[a]->id].natal_homing = this->natal_homing;
                females[areas_list[a]->id].area = areas_list[a];
                females[areas_list[a]->id].natal_area = this->natal_area;
                females[areas_list[a]->id].males = false;
                females[areas_list[a]->id].Initialize();
            }


        }

        void Prepare() {


            SN.resize(years * seasons * ages);
            SN_Biomass.resize(years * seasons * ages);
            C.resize(years * seasons * ages);
            C_Biomass.resize(years * seasons * ages);

            for (int a = 0; a < males.size(); a++) {
                males[areas_list[a]->id].Reset();
                females[areas_list[a]->id].Reset();
            }

            for (int i = 0; i < this->initial_numbers.size(); i++) {

                switch (this->initial_numbers[i].type) {
                    case MALE:
//                        std::cout << "Setting initial numbers for males in area " << this->initial_numbers[i].area_id << "\n";
                        this->males[this->initial_numbers[i].area_id].initial_numbers = this->initial_numbers[i].values;
                        break;

                    case FEMALE:
//                        std::cout << "Setting initial numbers for females in area " << this->initial_numbers[i].area_id << "\n";
                        this->females[this->initial_numbers[i].area_id].initial_numbers = this->initial_numbers[i].values;

                        break;

                }
            }



        }

        inline void IncrementTime(int& y, int& s) {
            if (s == this->seasons) {
                y += 1;
                s = 1;
            } else {
                s++;
            }
        }

        inline void DecrementTime(int& y, int& s) {
            if (s == 1) {
                y -= 1;
                s = seasons;
            } else {
                s--;
            }
        }

        inline void MoveFish(int year, int season) {

            int y = year;
            int s = season;
            //            IncrementTime(y, s);

            movement_model_iterator it = this->movement_models.find(year + 1);
            if (it != this->movement_models.end()) {



                int ss = season - 1;
                std::vector<std::vector<variable> >& male_probabilities = (*it).second->male_connectivity[ss];
                std::vector<std::vector<variable> >& female_probabilities = (*it).second->female_connectivity[ss];
                std::vector<std::vector<variable> >& rercruit_probabilities = (*it).second->recruit_connectivity[ss];
                //should be square
                for (int i = 0; i < male_probabilities.size(); i++) {
                    AreaPopulationInfo<REAL_T>& male_info_from = this->males[(i + 1)];
                    AreaPopulationInfo<REAL_T>& female_info_from = this->females[(i + 1)];
                    for (int j = 0; j < male_probabilities.size(); j++) {
                        AreaPopulationInfo<REAL_T>& male_info_to = this->males[(j + 1)];
                        AreaPopulationInfo<REAL_T>& female_info_to = this->females[(j + 1)];



                        if (i != j) {

                            variable tempm = rercruit_probabilities[i][j] * male_info_from.recruitment[year * this->seasons + (season - 1)];
                            variable tempf = rercruit_probabilities[i][j] * female_info_from.recruitment[year * this->seasons + (season - 1)];

                            male_info_to.redistributed_recruits[year * this->seasons + (season - 1)] += tempm;
                            female_info_to.redistributed_recruits[year * this->seasons + (season - 1)] += tempf;


                            for (int a = 0; a < this->ages; a++) {
                                male_info_from.emigrants[year * this->seasons * this->ages + (season - 1) * this->ages + a] +=
                                        male_probabilities[i][j] * male_info_from.N[year * this->seasons * this->ages + (season - 1) * this->ages + a];

                                male_info_to.imigrants[year * this->seasons * this->ages + (season - 1) * this->ages + a] +=
                                        male_probabilities[i][j] * male_info_from.N[year * this->seasons * this->ages + (season - 1) * this->ages + a];

                                female_info_from.emigrants[year * this->seasons * this->ages + (season - 1) * this->ages + a] +=
                                        female_probabilities[i][j] * female_info_from.N[year * this->seasons * this->ages + (season - 1) * this->ages + a];

                                female_info_to.imigrants[year * this->seasons * this->ages + (season - 1) * this->ages + a] +=
                                        female_probabilities[i][j] * female_info_from.N[year * this->seasons * this->ages + (season - 1) * this->ages + a];
                            }

                        }

                    }
                }



            } else {
                std::cout << "Configuration Error: Population " << this->id << " has no movement model defined for year " << (year + 1) << "\n";
                mas_log << "Configuration Error: Population " << this->id << " has no movement model defined for year " << (year + 1) << "\n";

            }

        }

        void InitializePopulationinAreas() {

            for (int a = 0; a < areas_list.size(); a++) {



                males[areas_list[a]->id].InitNumbers();
                females[areas_list[a]->id].InitNumbers();

            }

        }

        void Initialize() {

            for (int a = 0; a < areas_list.size(); a++) {

                males[areas_list[a]->id].Initialize();
                females[areas_list[a]->id].Initialize();

            }

        }

        void Show() {
            for (int a = 0; a < areas_list.size(); a++) {
                std::cout << males[areas_list[a]->id];
                std::cout << females[areas_list[a]->id];
            }


            std::cout << "Population: " << this->id << "\n";
            std::cout << "Total Catch Numbers At Age\n";
            for (int a = 0; a < this->ages; a++) {
                for (int y = 0; y < this->years; y++) {
                    for (int s = 1; s <= this->seasons; s++) {

                        std::cout << this->C[y * this->seasons * this->ages + (s - 1) * this->ages + a] << " ";
                    }
                }

                std::cout << "\n";
            }
            std::cout << "\n\n";

            std::cout << "Population: " << this->id << "\n";
            std::cout << "Total Catch Biomass At Age\n";
            for (int a = 0; a < this->ages; a++) {
                for (int y = 0; y < this->years; y++) {
                    for (int s = 1; s <= this->seasons; s++) {

                        std::cout << this->C_Biomass[y * this->seasons * this->ages + (s - 1) * this->ages + a] << " ";
                    }
                }

                std::cout << "\n";
            }
            std::cout << "\n\n";


            std::cout << "Population: " << this->id << "\n";
            std::cout << "Total Survey Numbers At Age\n";
            for (int a = 0; a < this->ages; a++) {
                for (int y = 0; y < this->years; y++) {
                    for (int s = 1; s <= this->seasons; s++) {

                        std::cout << this->SN[y * this->seasons * this->ages + (s - 1) * this->ages + a] << " ";
                    }
                }

                std::cout << "\n";
            }
            std::cout << "\n\n";

            std::cout << "Population: " << this->id << "\n";
            std::cout << "Total Survey Biomass At Age\n";
            for (int a = 0; a < this->ages; a++) {
                for (int y = 0; y < this->years; y++) {
                    for (int s = 1; s <= this->seasons; s++) {

                        std::cout << this->SN_Biomass[y * this->seasons * this->ages + (s - 1) * this->ages + a] << " ";
                    }

                }
                std::cout << "\n";
            }
            std::cout << "\n\n";
        }

        void Evaluate() {
            InitializePopulationinAreas();

            //
            //            if (this->move_fish_before_lh) {
            //                for (int y = 0; y < this->years; y++) {
            //                    for (int s = 1; s <= this->seasons; s++) {
            //
            //
            //                        for (int a = 0; a < areas_list.size(); a++) {
            //                            male_cohorts[areas_list[a]->id].Mortality(y, s);
            //                            female_cohorts[areas_list[a]->id].Mortality(y, s);
            //
            //
            //
            //                            male_cohorts[areas_list[a]->id].NumbersAtAge(y, s);
            //                            female_cohorts[areas_list[a]->id].NumbersAtAge(y, s);
            //
            //                            this->MoveFish(y, s);
            //
            //                        }
            //
            //
            //
            //
            //                        //                        for (int a = 0; a < areas_list.size(); a++) {
            //                        //                            male_cohorts[areas_list[a]->id].NumbersAtAge(y, s);
            //                        //                            female_cohorts[areas_list[a]->id].NumbersAtAge(y, s);
            //                        //                        }
            //
            //                        for (int a = 0; a < areas_list.size(); a++) {
            //                            male_cohorts[areas_list[a]->id].CatchAtAge(y, s);
            //                            female_cohorts[areas_list[a]->id].CatchAtAge(y, s);
            //                        }
            //
            //                    }
            //                }
            //            } else {
            for (int y = 0; y < this->years; y++) {
                for (int s = 1; s <= this->seasons; s++) {


                    for (int a = 0; a < areas_list.size(); a++) {




                        /******************************************
                         * Growth
                         *****************************************/
                        males[areas_list[a]->id].Growth(y, s);
                        females[areas_list[a]->id].Growth(y, s);


                        /******************************************
                         * Mortality
                         *****************************************/
                        males[areas_list[a]->id].Mortality(y, s);
                        females[areas_list[a]->id].Mortality(y, s);

                        /******************************************
                         * Numbers at Age
                         *****************************************/
                        males[areas_list[a]->id].NumbersAtAge(y, s);
                        females[areas_list[a]->id].NumbersAtAge(y, s);



                        /******************************************
                         * Spawning Biomass
                         *****************************************/
                        males[areas_list[a]->id].SpawningBiomass(y, s);
                        females[areas_list[a]->id].SpawningBiomass(y, s);

                        /******************************************
                         * Recruitment
                         *****************************************/
                        males[areas_list[a]->id].Recruitment(y, s);
                        females[areas_list[a]->id].Recruitment(y, s);

                        /******************************************
                         * Catch Numbers at Age
                         *****************************************/
                        males[areas_list[a]->id].CatchAtAge(y, s);
                        females[areas_list[a]->id].CatchAtAge(y, s);


                        /******************************************
                         * Survey Numbers at Age
                         *****************************************/
                        males[areas_list[a]->id].SurveyNumbersAtAge(y, s);
                        females[areas_list[a]->id].SurveyNumbersAtAge(y, s);


                    }
                    /******************************************
                     * Apply Movement
                     *****************************************/
                    this->MoveFish(y, s);
                }
            }


            for (int a = 0; a < areas_list.size(); a++) {
                /******************************************
                 * Push info to areas
                 *****************************************/
                males[areas_list[a]->id].PushToArea();
                females[areas_list[a]->id].PushToArea();
            }

            /**
             * Compute totals for this population
             */
            for (int y = 0; y < this->years; y++) {
                for (int s = 1; s <= this->seasons; s++) {


                    for (int al = 0; al < areas_list.size(); al++) {
                        for (int a = 0; a < this->ages; a++) {
                            this->C[y * this->seasons * this->ages + (s - 1) * this->ages + a] +=
                                    males[areas_list[al]->id].C[y * this->seasons * this->ages + (s - 1) * this->ages + a] +
                                    females[areas_list[al]->id].C[y * this->seasons * this->ages + (s - 1) * this->ages + a];

                            this->C_Biomass[y * this->seasons * this->ages + (s - 1) * this->ages + a] +=
                                    males[areas_list[al]->id].C_Biomass[y * this->seasons * this->ages + (s - 1) * this->ages + a] +
                                    females[areas_list[al]->id].C_Biomass[y * this->seasons * this->ages + (s - 1) * this->ages + a];

                            this->SN[y * this->seasons * this->ages + (s - 1) * this->ages + a] +=
                                    males[areas_list[al]->id].SN[y * this->seasons * this->ages + (s - 1) * this->ages + a] +
                                    females[areas_list[al]->id].SN[y * this->seasons * this->ages + (s - 1) * this->ages + a];

                            this->SN_Biomass[y * this->seasons * this->ages + (s - 1) * this->ages + a] +=
                                    males[areas_list[al]->id].SN_Biomass[y * this->seasons * this->ages + (s - 1) * this->ages + a] +
                                    females[areas_list[al]->id].SN_Biomass[y * this->seasons * this->ages + (s - 1) * this->ages + a];
                        }
                    }
                }
            }



        }

    };

    template<typename REAL_T>
    std::ostream& operator<<(std::ostream& out, const mas::Population<REAL_T>& pop) {
        out << "Population:\n";
        out << "Name: " << pop.name << "\n";
        out << "Id: " << pop.id << "\n";
        out << "Natal Area: " << pop.natal_area->id << "\n";
        out << "Movement Areas: ";
        for (int i = 0; i < pop.areas_list.size(); i++) {
            out << pop.areas_list[i]->id << " ";
        }
        out << "\n";

        return out;
    }

}

#endif /* POPULATION_HPP */

