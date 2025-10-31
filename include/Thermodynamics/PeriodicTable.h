/*
 *   This file is part of the OpenPhase (R) software library.
 *  
 *  Copyright (c) 2009-2025 Ruhr-Universitaet Bochum,
 *                Universitaetsstrasse 150, D-44801 Bochum, Germany
 *            AND 2018-2025 OpenPhase Solutions GmbH,
 *                Universitaetsstrasse 136, D-44799 Bochum, Germany.
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *     
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.

 *   File created :   2016
 *   Main contributors :   Matthias Stratmann
 *
 */

#ifndef PERIODICTABLE_H
#define PERIODICTABLE_H

#include "Includes.h"

namespace openphase
{

class PeriodicTable                                                             ///< Stores the properties of all chemical elements in the periodic table
{
 public:
    struct Element                                                              ///< Data structure for chemical element
    {
        std::string Name;
        std::string FullName;
        int AtomicNumber;
        double AtomicWeight;
        double MolarVolume;
    };

    PeriodicTable(void){};
    void Initialize(void)
    {
        Element temp;

        temp.Name = "H";
        temp.FullName = "Hydrogen";
        temp.AtomicNumber = 1;
        temp.AtomicWeight = 1.008;
        temp.MolarVolume  = 11.42E-6;

        Data.push_back(temp);

        temp.Name = "HE";
        temp.FullName = "Helium";
        temp.AtomicNumber = 2;
        temp.AtomicWeight = 4.003;
        temp.MolarVolume  = 21.00E-6;

        Data.push_back(temp);

        temp.Name = "LI";
        temp.FullName = "Lithium";
        temp.AtomicNumber = 3;
        temp.AtomicWeight = 6.940;
        temp.MolarVolume  = 13.02E-6;

        Data.push_back(temp);

        temp.Name = "BE";
        temp.FullName = "Beryllium";
        temp.AtomicNumber = 4;
        temp.AtomicWeight = 9.012;
        temp.MolarVolume  = 4.85E-6;

        Data.push_back(temp);

        temp.Name = "B";
        temp.FullName = "Boron";
        temp.AtomicNumber = 5;
        temp.AtomicWeight = 10.810;
        temp.MolarVolume  = 4.39E-6;

        Data.push_back(temp);

        temp.Name = "C";
        temp.FullName = "Carbon";
        temp.AtomicNumber = 6;
        temp.AtomicWeight = 12.011;
        temp.MolarVolume  = 5.31E-6;

        Data.push_back(temp);

        temp.Name = "N";
        temp.FullName = "Nitrogen";
        temp.AtomicNumber = 7;
        temp.AtomicWeight = 14.007;
        temp.MolarVolume  = 13.54E-6;

        Data.push_back(temp);

        temp.Name = "O";
        temp.FullName = "Oxygen";
        temp.AtomicNumber = 8;
        temp.AtomicWeight = 15.999;
        temp.MolarVolume  = 17.36E-6;

        Data.push_back(temp);

        temp.Name = "F";
        temp.FullName = "Fluorine";
        temp.AtomicNumber = 9;
        temp.AtomicWeight = 18.998;
        temp.MolarVolume  = 11.20E-6;

        Data.push_back(temp);

        temp.Name = "NE";
        temp.FullName = "Neon";
        temp.AtomicNumber = 10;
        temp.AtomicWeight = 20.180;
        temp.MolarVolume  = 13.23E-6;

        Data.push_back(temp);

        temp.Name = "NA";
        temp.FullName = "Sodium";
        temp.AtomicNumber = 11;
        temp.AtomicWeight = 22.990;
        temp.MolarVolume  = 23.78E-6;

        Data.push_back(temp);

        temp.Name = "MG";
        temp.FullName = "Magnesium";
        temp.AtomicNumber = 12;
        temp.AtomicWeight = 24.305;
        temp.MolarVolume  = 14.00E-6;

        Data.push_back(temp);

        temp.Name = "AL";
        temp.FullName = "Aluminium";
        temp.AtomicNumber = 13;
        temp.AtomicWeight = 26.982;
        temp.MolarVolume  = 10.00E-6;

        Data.push_back(temp);

        temp.Name = "SI";
        temp.FullName = "Silicon";
        temp.AtomicNumber = 14;
        temp.AtomicWeight = 28.085;
        temp.MolarVolume  = 12.06E-6;

        Data.push_back(temp);

        temp.Name = "P";
        temp.FullName = "Phosphorus";
        temp.AtomicNumber = 15;
        temp.AtomicWeight = 30.974;
        temp.MolarVolume  = 17.02E-6;

        Data.push_back(temp);

        temp.Name = "S";
        temp.FullName = "Sulfur";
        temp.AtomicNumber = 16;
        temp.AtomicWeight = 32.060;
        temp.MolarVolume  = 15.53E-6;

        Data.push_back(temp);

        temp.Name = "CL";
        temp.FullName = "Chlorine";
        temp.AtomicNumber = 17;
        temp.AtomicWeight = 35.450;
        temp.MolarVolume  = 17.39E-6;

        Data.push_back(temp);

        temp.Name = "AR";
        temp.FullName = "Argon";
        temp.AtomicNumber = 18;
        temp.AtomicWeight = 39.948;
        temp.MolarVolume  = 22.56E-6;

        Data.push_back(temp);

        temp.Name = "K";
        temp.FullName = "Potassium";
        temp.AtomicNumber = 19;
        temp.AtomicWeight = 39.098;
        temp.MolarVolume  = 45.94E-6;

        Data.push_back(temp);

        temp.Name = "CA";
        temp.FullName = "Calcium";
        temp.AtomicNumber = 20;
        temp.AtomicWeight = 40.078;
        temp.MolarVolume  = 26.20E-6;

        Data.push_back(temp);

        temp.Name = "SC";
        temp.FullName = "Scandium";
        temp.AtomicNumber = 21;
        temp.AtomicWeight = 44.956;
        temp.MolarVolume  = 15.00E-6;

        Data.push_back(temp);

        temp.Name = "TI";
        temp.FullName = "Titanium";
        temp.AtomicNumber = 22;
        temp.AtomicWeight = 47.867;
        temp.MolarVolume  = 10.64E-6;

        Data.push_back(temp);

        temp.Name = "V";
        temp.FullName = "Vanadium";
        temp.AtomicNumber = 23;
        temp.AtomicWeight = 50.942;
        temp.MolarVolume  = 8.32E-6;

        Data.push_back(temp);

        temp.Name = "CR";
        temp.FullName = "Chromium";
        temp.AtomicNumber = 24;
        temp.AtomicWeight = 51.996;
        temp.MolarVolume  = 7.23E-6;

        Data.push_back(temp);

        temp.Name = "MN";
        temp.FullName = "Manganese";
        temp.AtomicNumber = 25;
        temp.AtomicWeight = 54.938;
        temp.MolarVolume  = 7.35E-6;

        Data.push_back(temp);

        temp.Name = "FE";
        temp.FullName = "Iron";
        temp.AtomicNumber = 26;
        temp.AtomicWeight = 55.845;
        temp.MolarVolume  = 7.09E-6;

        Data.push_back(temp);

        temp.Name = "CO";
        temp.FullName = "Cobalt";
        temp.AtomicNumber = 27;
        temp.AtomicWeight = 58.933;
        temp.MolarVolume  = 6.67E-6;

        Data.push_back(temp);

        temp.Name = "NI";
        temp.FullName = "Nickel";
        temp.AtomicNumber = 28;
        temp.AtomicWeight = 58.693;
        temp.MolarVolume  = 6.59E-6;

        Data.push_back(temp);

        temp.Name = "CU";
        temp.FullName = "Copper";
        temp.AtomicNumber = 29;
        temp.AtomicWeight = 63.546;
        temp.MolarVolume  = 7.11E-6;

        Data.push_back(temp);

        temp.Name = "ZN";
        temp.FullName = "Zinc";
        temp.AtomicNumber = 30;
        temp.AtomicWeight = 65.380;
        temp.MolarVolume  = 9.16E-6;

        Data.push_back(temp);

        temp.Name = "GA";
        temp.FullName = "Gallium";
        temp.AtomicNumber = 31;
        temp.AtomicWeight = 69.723;
        temp.MolarVolume  = 11.80E-6;

        Data.push_back(temp);

        temp.Name = "GE";
        temp.FullName = "Germanium";
        temp.AtomicNumber = 32;
        temp.AtomicWeight = 72.630;
        temp.MolarVolume  = 13.63E-6;

        Data.push_back(temp);

        temp.Name = "AS";
        temp.FullName = "Arsenic";
        temp.AtomicNumber = 33;
        temp.AtomicWeight = 74.922;
        temp.MolarVolume  = 12.95E-6;

        Data.push_back(temp);

        temp.Name = "SE";
        temp.FullName = "Selenium";
        temp.AtomicNumber = 34;
        temp.AtomicWeight = 78.971;
        temp.MolarVolume  = 16.42E-6;

        Data.push_back(temp);

        temp.Name = "BR";
        temp.FullName = "Bromine";
        temp.AtomicNumber = 35;
        temp.AtomicWeight = 79.904;
        temp.MolarVolume  = 19.78E-6;

        Data.push_back(temp);

        temp.Name = "KR";
        temp.FullName = "Krypton";
        temp.AtomicNumber = 36;
        temp.AtomicWeight = 83.798;
        temp.MolarVolume  = 27.99E-6;

        Data.push_back(temp);

        temp.Name = "RB";
        temp.FullName = "Rubidium";
        temp.AtomicNumber = 37;
        temp.AtomicWeight = 85.468;
        temp.MolarVolume  = 55.76E-6;

        Data.push_back(temp);

        temp.Name = "SR";
        temp.FullName = "Strontium";
        temp.AtomicNumber = 38;
        temp.AtomicWeight = 87.620;
        temp.MolarVolume  = 33.94E-6;

        Data.push_back(temp);

        temp.Name = "Y";
        temp.FullName = "Yttrium";
        temp.AtomicNumber = 39;
        temp.AtomicWeight = 88.906;
        temp.MolarVolume  = 19.88E-6;

        Data.push_back(temp);

        temp.Name = "ZR";
        temp.FullName = "Zirconium";
        temp.AtomicNumber = 40;
        temp.AtomicWeight = 91.224;
        temp.MolarVolume  = 14.02E-6;

        Data.push_back(temp);

        temp.Name = "NB";
        temp.FullName = "Niobium";
        temp.AtomicNumber = 41;
        temp.AtomicWeight = 92.906;
        temp.MolarVolume  = 10.83E-6;

        Data.push_back(temp);

        temp.Name = "MO";
        temp.FullName = "Molybdenum";
        temp.AtomicNumber = 42;
        temp.AtomicWeight = 95.951;
        temp.MolarVolume  = 9.38E-6;

        Data.push_back(temp);

        temp.Name = "TC";
        temp.FullName = "Technetium";
        temp.AtomicNumber = 43;
        temp.AtomicWeight = 98.906;
        temp.MolarVolume  = 8.63E-6;

        Data.push_back(temp);

        temp.Name = "RU";
        temp.FullName = "Ruthenium";
        temp.AtomicNumber = 44;
        temp.AtomicWeight = 101.072;
        temp.MolarVolume  = 8.17E-6;

        Data.push_back(temp);

        temp.Name = "RH";
        temp.FullName = "Rhodium";
        temp.AtomicNumber = 45;
        temp.AtomicWeight = 102.906;
        temp.MolarVolume  = 8.28E-6;

        Data.push_back(temp);

        temp.Name = "PD";
        temp.FullName = "Palladium";
        temp.AtomicNumber = 46;
        temp.AtomicWeight = 106.421;
        temp.MolarVolume  = 8.56E-6;

        Data.push_back(temp);

        temp.Name = "AG";
        temp.FullName = "Silver";
        temp.AtomicNumber = 47;
        temp.AtomicWeight = 107.868;
        temp.MolarVolume  = 10.27E-6;

        Data.push_back(temp);

        temp.Name = "CD";
        temp.FullName = "Cadmium";
        temp.AtomicNumber = 48;
        temp.AtomicWeight = 112.414;
        temp.MolarVolume  = 13.00E-6;

        Data.push_back(temp);

        temp.Name = "IN";
        temp.FullName = "Indium";
        temp.AtomicNumber = 49;
        temp.AtomicWeight = 114.818;
        temp.MolarVolume  = 15.76E-6;

        Data.push_back(temp);

        temp.Name = "SN";
        temp.FullName = "Tin";
        temp.AtomicNumber = 50;
        temp.AtomicWeight = 118.710;
        temp.MolarVolume  = 16.29E-6;

        Data.push_back(temp);

        temp.Name = "SB";
        temp.FullName = "Antimony";
        temp.AtomicNumber = 51;
        temp.AtomicWeight = 121.760;
        temp.MolarVolume  = 18.19E-6;

        Data.push_back(temp);

        temp.Name = "TE";
        temp.FullName = "Tellurium";
        temp.AtomicNumber = 52;
        temp.AtomicWeight = 127.603;
        temp.MolarVolume  = 20.46E-6;

        Data.push_back(temp);

        temp.Name = "I";
        temp.FullName = "Iodine";
        temp.AtomicNumber = 53;
        temp.AtomicWeight = 126.904;
        temp.MolarVolume  = 25.72E-6;

        Data.push_back(temp);

        temp.Name = "XE";
        temp.FullName = "Xenon";
        temp.AtomicNumber = 54;
        temp.AtomicWeight = 131.293;
        temp.MolarVolume  = 35.92E-6;

        Data.push_back(temp);

        temp.Name = "CS";
        temp.FullName = "Caesium";
        temp.AtomicNumber = 55;
        temp.AtomicWeight = 132.905;
        temp.MolarVolume  = 70.94E-6;

        Data.push_back(temp);

        temp.Name = "BA";
        temp.FullName = "Barium";
        temp.AtomicNumber = 56;
        temp.AtomicWeight = 137.327;
        temp.MolarVolume  = 38.16E-6;

        Data.push_back(temp);

        temp.Name = "LA";
        temp.FullName = "Lanthanum";
        temp.AtomicNumber = 57;
        temp.AtomicWeight = 138.905;
        temp.MolarVolume  = 22.39E-6;

        Data.push_back(temp);

        temp.Name = "CE";
        temp.FullName = "Cerium";
        temp.AtomicNumber = 58;
        temp.AtomicWeight = 140.116;
        temp.MolarVolume  = 20.69E-6;

        Data.push_back(temp);

        temp.Name = "PR";
        temp.FullName = "Praseodymium";
        temp.AtomicNumber = 59;
        temp.AtomicWeight = 140.908;
        temp.MolarVolume  = 20.80E-6;

        Data.push_back(temp);

        temp.Name = "ND";
        temp.FullName = "Neodymium";
        temp.AtomicNumber = 60;
        temp.AtomicWeight = 144.242;
        temp.MolarVolume  = 20.59E-6;

        Data.push_back(temp);

        temp.Name = "PM";
        temp.FullName = "Promethium";
        temp.AtomicNumber = 61;
        temp.AtomicWeight = 146.915;
        temp.MolarVolume  = 20.10E-6;

        Data.push_back(temp);

        temp.Name = "SM";
        temp.FullName = "Samarium";
        temp.AtomicNumber = 62;
        temp.AtomicWeight = 150.362;
        temp.MolarVolume  = 19.98E-6;

        Data.push_back(temp);

        temp.Name = "EU";
        temp.FullName = "Europium";
        temp.AtomicNumber = 63;
        temp.AtomicWeight = 151.964;
        temp.MolarVolume  = 28.97E-6;

        Data.push_back(temp);

        temp.Name = "GD";
        temp.FullName = "Gadolinium";
        temp.AtomicNumber = 64;
        temp.AtomicWeight = 157.253;
        temp.MolarVolume  = 19.90E-6;

        Data.push_back(temp);

        temp.Name = "TB";
        temp.FullName = "Terbium";
        temp.AtomicNumber = 65;
        temp.AtomicWeight = 158.925;
        temp.MolarVolume  = 19.30E-6;

        Data.push_back(temp);

        temp.Name = "DY";
        temp.FullName = "Dysprosium";
        temp.AtomicNumber = 66;
        temp.AtomicWeight = 162.500;
        temp.MolarVolume  = 19.01E-6;

        Data.push_back(temp);

        temp.Name = "HO";
        temp.FullName = "Holmium";
        temp.AtomicNumber = 67;
        temp.AtomicWeight = 164.930;
        temp.MolarVolume  = 18.74E-6;

        Data.push_back(temp);

        temp.Name = "ER";
        temp.FullName = "Erbium";
        temp.AtomicNumber = 68;
        temp.AtomicWeight = 167.259;
        temp.MolarVolume  = 18.46E-6;

        Data.push_back(temp);

        temp.Name = "TM";
        temp.FullName = "Thulium";
        temp.AtomicNumber = 69;
        temp.AtomicWeight = 168.934;
        temp.MolarVolume  = 19.10E-6;

        Data.push_back(temp);

        temp.Name = "YB";
        temp.FullName = "Ytterbium";
        temp.AtomicNumber = 70;
        temp.AtomicWeight = 173.045;
        temp.MolarVolume  = 24.84E-6;

        Data.push_back(temp);

        temp.Name = "LU";
        temp.FullName = "Lutetium";
        temp.AtomicNumber = 71;
        temp.AtomicWeight = 174.967;
        temp.MolarVolume  = 17.78E-6;

        Data.push_back(temp);

        temp.Name = "HF";
        temp.FullName = "Hafnium";
        temp.AtomicNumber = 72;
        temp.AtomicWeight = 178.492;
        temp.MolarVolume  = 13.44E-6;

        Data.push_back(temp);

        temp.Name = "TA";
        temp.FullName = "Tantalum";
        temp.AtomicNumber = 73;
        temp.AtomicWeight = 180.948;
        temp.MolarVolume  = 10.85E-6;

        Data.push_back(temp);

        temp.Name = "W";
        temp.FullName = "Tungsten";
        temp.AtomicNumber = 74;
        temp.AtomicWeight = 183.841;
        temp.MolarVolume  = 9.47E-6;

        Data.push_back(temp);

        temp.Name = "RE";
        temp.FullName = "Rhenium";
        temp.AtomicNumber = 75;
        temp.AtomicWeight = 186.207;
        temp.MolarVolume  = 8.86E-6;

        Data.push_back(temp);

        temp.Name = "OS";
        temp.FullName = "Osmium";
        temp.AtomicNumber = 76;
        temp.AtomicWeight = 190.233;
        temp.MolarVolume  = 8.42E-6;

        Data.push_back(temp);

        temp.Name = "IR";
        temp.FullName = "Iridium";
        temp.AtomicNumber = 77;
        temp.AtomicWeight = 192.217;
        temp.MolarVolume  = 8.52E-6;

        Data.push_back(temp);

        temp.Name = "PT";
        temp.FullName = "Platinum";
        temp.AtomicNumber = 78;
        temp.AtomicWeight = 195.084;
        temp.MolarVolume  = 9.09E-6;

        Data.push_back(temp);

        temp.Name = "AU";
        temp.FullName = "Gold";
        temp.AtomicNumber = 79;
        temp.AtomicWeight = 196.967;
        temp.MolarVolume  = 10.21E-6;

        Data.push_back(temp);

        temp.Name = "HG";
        temp.FullName = "Mercury";
        temp.AtomicNumber = 80;
        temp.AtomicWeight = 200.592;
        temp.MolarVolume  = 14.09E-6;

        Data.push_back(temp);

        temp.Name = "TL";
        temp.FullName = "Thallium";
        temp.AtomicNumber = 81;
        temp.AtomicWeight = 204.382;
        temp.MolarVolume  = 17.22E-6;

        Data.push_back(temp);

        temp.Name = "PB";
        temp.FullName = "Lead";
        temp.AtomicNumber = 82;
        temp.AtomicWeight = 207.210;
        temp.MolarVolume  = 18.26E-6;

        Data.push_back(temp);

        temp.Name = "BI";
        temp.FullName = "Bismuth";
        temp.AtomicNumber = 83;
        temp.AtomicWeight = 208.980;
        temp.MolarVolume  = 21.31E-6;

        Data.push_back(temp);

        temp.Name = "PO";
        temp.FullName = "Polonium";
        temp.AtomicNumber = 84;
        temp.AtomicWeight = 209.980;
        temp.MolarVolume  = 22.97E-6;

        Data.push_back(temp);

        temp.Name = "AT";
        temp.FullName = "Astatine";
        temp.AtomicNumber = 85;
        temp.AtomicWeight = 209.987;
        temp.MolarVolume  = -1.0;

        Data.push_back(temp);

        temp.Name = "RN";
        temp.FullName = "Radon";
        temp.AtomicNumber = 86;
        temp.AtomicWeight = 222.000;
        temp.MolarVolume  = 50.50E-6;

        Data.push_back(temp);

        temp.Name = "FR";
        temp.FullName = "Francium";
        temp.AtomicNumber = 87;
        temp.AtomicWeight = 223.020;
        temp.MolarVolume  = -1.0;

        Data.push_back(temp);

        temp.Name = "RA";
        temp.FullName = "Radium";
        temp.AtomicNumber = 88;
        temp.AtomicWeight = 226.025;
        temp.MolarVolume  = 41.09E-6;

        Data.push_back(temp);

        temp.Name = "AC";
        temp.FullName = "Actinium";
        temp.AtomicNumber = 89;
        temp.AtomicWeight = 227.028;
        temp.MolarVolume  = 22.55E-6;

        Data.push_back(temp);

        temp.Name = "TH";
        temp.FullName = "Thorium";
        temp.AtomicNumber = 90;
        temp.AtomicWeight = 232.038;
        temp.MolarVolume  = 19.80E-6;

        Data.push_back(temp);

        temp.Name = "PA";
        temp.FullName = "Protactinium";
        temp.AtomicNumber = 91;
        temp.AtomicWeight = 231.036;
        temp.MolarVolume  = 15.18E-6;

        Data.push_back(temp);

        temp.Name = "U";
        temp.FullName = "Uranium";
        temp.AtomicNumber = 92;
        temp.AtomicWeight = 238.029;
        temp.MolarVolume  = 12.49E-6;

        Data.push_back(temp);

        temp.Name = "NP";
        temp.FullName = "Neptunium";
        temp.AtomicNumber = 93;
        temp.AtomicWeight = 237.048;
        temp.MolarVolume  = 11.59E-6;

        Data.push_back(temp);

        temp.Name = "PU";
        temp.FullName = "Plutonium";
        temp.AtomicNumber = 94;
        temp.AtomicWeight = 244.064;
        temp.MolarVolume  = 12.29E-6;

        Data.push_back(temp);

        temp.Name = "AM";
        temp.FullName = "Americium";
        temp.AtomicNumber = 95;
        temp.AtomicWeight = 243.061;
        temp.MolarVolume  = 17.78E-6;

        Data.push_back(temp);

        temp.Name = "CM";
        temp.FullName = "Curium";
        temp.AtomicNumber = 96;
        temp.AtomicWeight = 247.070;
        temp.MolarVolume  = 18.05E-6;

        Data.push_back(temp);

        temp.Name = "BK";
        temp.FullName = "Berkelium";
        temp.AtomicNumber = 97;
        temp.AtomicWeight = 247.0;
        temp.MolarVolume  = 16.84E-6;

        Data.push_back(temp);

        temp.Name = "CF";
        temp.FullName = "Californium";
        temp.AtomicNumber = 98;
        temp.AtomicWeight = 251.0;
        temp.MolarVolume  = 16.50E-6;

        Data.push_back(temp);

        temp.Name = "ES";
        temp.FullName = "Einsteinium";
        temp.AtomicNumber = 99;
        temp.AtomicWeight = 252.0;
        temp.MolarVolume  = -1.0;

        Data.push_back(temp);

        temp.Name = "FM";
        temp.FullName = "Fermium";
        temp.AtomicNumber = 100;
        temp.AtomicWeight = 257.095;
        temp.MolarVolume  = -1.0;

        Data.push_back(temp);

        temp.Name = "MD";
        temp.FullName = "Mendelevium";
        temp.AtomicNumber = 101;
        temp.AtomicWeight = 258.0;
        temp.MolarVolume  = -1.0;

        Data.push_back(temp);

        temp.Name = "NO";
        temp.FullName = "Nobelium";
        temp.AtomicNumber = 102;
        temp.AtomicWeight = 259.0;
        temp.MolarVolume  = -1.0;

        Data.push_back(temp);

        temp.Name = "LR";
        temp.FullName = "Lawrencium";
        temp.AtomicNumber = 103;
        temp.AtomicWeight = 266.0;
        temp.MolarVolume  = -1.0;

        Data.push_back(temp);

        temp.Name = "RD";
        temp.FullName = "Rutherfordium";
        temp.AtomicNumber = 104;
        temp.AtomicWeight = 261.109;
        temp.MolarVolume  = -1.0;

        Data.push_back(temp);

        temp.Name = "DB";
        temp.FullName = "Dubnium";
        temp.AtomicNumber = 105;
        temp.AtomicWeight = 262.114;
        temp.MolarVolume  = -1.0;

        Data.push_back(temp);

        temp.Name = "SG";
        temp.FullName = "Seaborgium";
        temp.AtomicNumber = 106;
        temp.AtomicWeight = 263.118;
        temp.MolarVolume  = -1.0;

        Data.push_back(temp);

        temp.Name = "BH";
        temp.FullName = "Bohrium";
        temp.AtomicNumber = 107;
        temp.AtomicWeight = 262.123;
        temp.MolarVolume  = -1.0;

        Data.push_back(temp);

        temp.Name = "HS";
        temp.FullName = "Hassium";
        temp.AtomicNumber = 108;
        temp.AtomicWeight = 265.0;
        temp.MolarVolume  = -1.0;

        Data.push_back(temp);

        temp.Name = "MT";
        temp.FullName = "Meitnerium";
        temp.AtomicNumber = 109;
        temp.AtomicWeight = 268.0;
        temp.MolarVolume  = -1.0;

        Data.push_back(temp);

        temp.Name = "DS";
        temp.FullName = "Darmstadtium";
        temp.AtomicNumber = 110;
        temp.AtomicWeight = 281.;
        temp.MolarVolume  = -1.;

        Data.push_back(temp);

        temp.Name = "RG";
        temp.FullName = "Roentgenium";
        temp.AtomicNumber = 111;
        temp.AtomicWeight = 280.0;
        temp.MolarVolume  = -1.0;

        Data.push_back(temp);

        temp.Name = "CN";
        temp.FullName = "Copernicium";
        temp.AtomicNumber = 112;
        temp.AtomicWeight = 277.0;
        temp.MolarVolume  = -1.0;

        Data.push_back(temp);

        temp.Name = "NH";
        temp.FullName = "Nihonium";
        temp.AtomicNumber = 113;
        temp.AtomicWeight = 287.;
        temp.MolarVolume  = -1.;

        Data.push_back(temp);

        temp.Name = "FL";
        temp.FullName = "Flerovium";
        temp.AtomicNumber = 114;
        temp.AtomicWeight = 289.0;
        temp.MolarVolume  = -1.0;

        Data.push_back(temp);

        temp.Name = "MC";
        temp.FullName = "Moscovium";
        temp.AtomicNumber = 115;
        temp.AtomicWeight = 288.0;
        temp.MolarVolume  = -1.0;

        Data.push_back(temp);

        temp.Name = "LV";
        temp.FullName = "Livermorium";
        temp.AtomicNumber = 116;
        temp.AtomicWeight = 293.0;
        temp.MolarVolume  = -1.0;

        Data.push_back(temp);

        temp.Name = "TS";
        temp.FullName = "Tennessine";
        temp.AtomicNumber = 117;
        temp.AtomicWeight = 292.0;
        temp.MolarVolume  = -1.0;

        Data.push_back(temp);

        temp.Name = "OG";
        temp.FullName = "Oganesson";
        temp.AtomicNumber = 118;
        temp.AtomicWeight = 294.0;
        temp.MolarVolume  = -1.0;

        Data.push_back(temp);
    }

    Element GetData(std::string Name) const
    {
        std::transform(Name.begin(), Name.end(), Name.begin(), ::toupper);      //Change "Element" to upper cases only
        for(size_t n = 0; n < Data.size(); n++)
        {
            if(Name == Data[n].Name)
            {
                return Data[n];
            }
        }
        Element temp;
        temp.Name = Name;
        temp.FullName = Name;
        temp.AtomicNumber = 100;
        temp.AtomicWeight = 100.0;
        temp.MolarVolume  = 1E-5;
        //std::cout << "Could not find Element " << Name
        //<< " in PeriodicTable.h! Used Custom Element instead!" << std::endl;
        return temp;
    }

 protected:
 private:
    std::vector<Element> Data;
};
}
#endif//PERIODICTABLE_H
