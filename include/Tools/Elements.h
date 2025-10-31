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

 *   File created :   2015
 *   Main contributors :   Rajendran Mohan Kumar
 *
 */

#ifndef ELEMENTS_H
#define ELEMENTS_H

#define ANGSTROM "&#x212B;"     // Unicode for angstrom sign
#define DEGREE "&#176;"         // Unicode for degree sign

#include "Includes.h"

namespace openphase {

class ConsoleOutput;
class Elements
{
public:

    const static std::string thisclassname = "Elements";

    const static std::vector<std::string> headers =
    { "Name","Symbol","Atomic Number","Atomic Weight (u.m)","Density (g/cm<sup>3</sup>) ","Melting Point (K)","Boiling Point (K)","Atomic Radius (pm)","Covalent Radius (pm)","Ionic Radius (pm)","Atomic Volume (cm<sup>3</sup>/mol)","Specific Heat (@20""DEGREE""C J/g mol)","Fusion Heat (kJ/mol)","Evaporation Heat (kJ/mol)","Termal Conductivity (@25""DEGREE""C W/m K)","Debye temperature (K)","Pauling Negativity Number","First Ionizing Energy (kJ/mol)","Oxidation states","Electronic Configuration","Lattice structure","Lattice constant (""ANGSTROM"")","Lattice c/a ratio","Appearance","Discovery date","Discovered by","Named after","X","Y" };

    const static std::vector<std::vector<std::string>> periodicTable = {
    { "Hydrogen","H","1","1.00794","0.0708 (@ -253°C)","14.01","20.28","79","32","154 (-1e)","14.1","14.267 (H-H)","0.117 (H-H)","0.904 (H-H)","0.1815","110","2.2","1311.3","1, -1","1s<sup>1</sup>","HEX","3.75","1.731","Colorless, odorless, tasteless gas","1766 (England)","Henry Cavendish","Greek: hydro (water) and genes (generate)","1","1" },
    { "Helium","He","2","4.002602","0.147 (@ -270°C)","0.95","4.216","0","n/a","93","31.8","5.188","n/a","0.08","0.152","n/a","n/a","2361.3","n/a","1s<sup>2</sup>","HEX","3.57","1.633","Inert, colorless, odorless, tasteless gas","1895 (Scotland/Sweden)","Sir William Ramsey, Nils Langet, P.T.Cleve","Greek: helios (sun).","18","1" },
    { "Lithium","Li","3","6.941","0.534","553.69","1118.15","155","163","68 (+1e)","13.1","3.489","2.89","148","84.8","400","0.98","519.9","1","[He] 2s<sup>1</sup>","BCC","3.49","n/a","Soft, silvery-white metal","1817 (Sweden)","Johann Arfwedson","Greek: lithos (stone).","1","2" },
    { "Beryllium","Be","4","9.01218","1.848","1551","3243","112","90","35 (+2e)","5","1.824","12.21","309","201","1000","1.57","898.8","2","[He] 2s<sup>2</sup>","HEX","2.29","1.567","Hard, brittle, steel-gray metal","1798 (France)","Louis-Nicholas Vauquelin","Greek: beryllos, 'beryl' (a mineral).","2","2" },
    { "Boron","B","5","10.811","2.34","2573","3931","98","82","23 (+3e)","4.6","1.025","23.6","504.5","27.4","1250","2.04","800.2","3","[He] 2s<sup>2</sup> 2p<sup>1</sup>","TET","8.73","0.576","Hard, brittle, lustrous black semimetal","1808 (England/France)","Sir H. Davy, J.L. Gay-Lussac, L.J. Thenard","The Arabic and Persian words for borax.","13","2" },
    { "Carbon","C","6","12.011","2.25 (graphite)","3820","5100","91","77","16 (+4e) 260 (-4e)","5.3","0.711","n/a","n/a","1.59","1860","2.55","1085.7","4, 2, -4","[He] 2s<sup>2</sup> 2p<sup>2</sup>","DIA","3.57","n/a","Dense, Black","n/a (Unknown)","Known to the ancients","Latin: carbo, (charcoal).","14","2" },
    { "Nitrogen","N","7","14.00674","0.808 (@ -195.8&#176;C)","63.29","77.4","92","75","13 (+5e) 171 (-3e)","17.3","1.042 (N-N)","n/a","n/a","0.026","n/a","3.04","1401.5","5, 4, 3, 2, -3","[He] 2s<sup>2</sup> 2p<sup>3</sup>","HEX","4.039","1.651","Colorless, odorless, tasteless, and generally inert gas","1772 (Scotland)","Daniel Rutherford","Greek: nitron and genes, (soda forming).","15","2" },
    { "Oxygen","O","8","15.9994","1.149 (@ -183&#176;C)","54.8","90.19","n/a","73","132 (-2e)","14","0.916 (O-O)","n/a","n/a","0.027","n/a","3.44","1313.1","-2, -1","[He] 2s<sup>2</sup> 2p<sup>4</sup>","CUB","6.83","n/a","Colorless, odorless, tasteless gas; pale blue liquid","1774 (England/Sweden)","Joseph Priestly, Carl Wilhelm Scheele","Greek: oxys and genes, (acid former).","16","2" },
    { "Fluorine","F","9","18.998403","1.108 (@ -189&#176;C)","53.53","85.01","n/a","72","133 (-1e)","17.1","0.824 (F-F)","0.51 (F-F)","6.54 (F-F)","0.028","n/a","3.98","1680","-1","[He] 2s<sup>2</sup> 2p<sup>5</sup>","MCL","n/a","n/a","Greenish-yellow, pungent, corrosive gas","1886 (France)","Henri Moissan","Latin: fluere (flow).","17","2" },
    { "Neon","Ne","10","20.1797","1.204 (@ -246&#176;C)","48","27.1","n/a","71","n/a","16.8","1.029","n/a","1.74","-0.0493","63","0","2079.4","n/a","[He] 2s<sup>2</sup> 2p<sup>6</sup>","FCC","4.43","n/a","Colorless, odorless, tasteless gas","1898 (England)","Sir William Ramsey, M.W. Travers","Greek: neos (new).","18","2" },
    { "Sodium","Na","11","22.989768","0.971","370.96","1156.1","190","154","97 (+1e)","23.7","1.222","2.64","97.9","142","150","0.93","495.6","1","[Ne] 3s<sup>1</sup>","BCC","4.23","n/a","Soft, silvery-white metal","1807 (England)","Sir Humphrey Davy","Medieval Latin: sodanum, (headache remedy); symbol from Latin natrium, (sodium carbonate).","1","3" },
    { "Magnesium","Mg","12","24.305","1.738","922","1363","160","136","66 (+2e)","14","1.025","9.2","131.8","156","318","1.31","737.3","2","[Ne] 3s<sup>2</sup>","HEX","3.21","1.624","Lightweight, malleable, silvery-white metal","1808 (England)","Sir Humphrey Davy","Magnesia, ancient city in district of Thessaly, Greece.","2","3" },
    { "Aluminum","Al","13","26.981539","2.6989","933.5","2740","143","118","51 (+3e)","10","0.9","10.75","284.1","237","394","1.61","577.2","3","[Ne] 3s<sup>2</sup> 3p<sup>1</sup>","FCC","4.05","n/a","Soft, lightweight, silvery-white metal","1825 (Denmark)","Hans Christian Oersted","Latin: alumen, aluminis, (alum).","13","3" },
    { "Silicon","Si","14","28.0855","2.33","1683","2628","132","111","42 (+4e) 271  (-4e)","12.1","0.703","50.6","383","149","625","1.9","786","4, -4","[Ne] 3s<sup>2</sup> 3p<sup>2</sup>","DIA","5.43","n/a","Amorphous form is brown powder; crystalline form has a gray","1824 (Sweden)","Jons Jacob Berzelius","Latin: silex, silicus, (flint).","14","3" },
    { "Phosphorus","P","15","30.973762","1.82 (white phosphorus)","317.3","553","128","106","35 (+5e) 212 (-3e)","17","0.757","2.51","49.8","-0.236","n/a","2.19","1011.2","5, 3, -3","[Ne] 3s<sup>2</sup> 3p<sup>3</sup>","CUB","7.17","n/a","The most common white form is a waxy, phosphorescent solid","1669 (Germany)","Hennig Brand","Greek: phosphoros, (bringer of light).","15","3" },
    { "Sulfur","S","16","32.066","2.07","386","717.824","127","102","30 (+6e) 184 (-2e)","15.5","0.732","1.23","10.5","0.27","n/a","2.58","999","6, 4, 2, -2","[Ne] 3s<sup>2</sup> 3p<sup>4</sup>","ORC","10.47","n/a","Tasteless, odorless, light-yellow, brittle solid","n/a (Unknown)","Known to the ancients.","Latin: sulphur (brimstone).","16","3" },
    { "Chlorine","Cl","17","35.4527","1.56 (@ -33.6Â°C)","172.2","238.6","n/a","99","27 (+7e) 181 (-1e)","18.7","0.477 (Cl-Cl)","6.41 (Cl-Cl)","20.41 (Cl-Cl)","0.009","n/a","3.16","1254.9","7, 5, 3, 1, -1","[Ne] 3s<sup>2</sup> 3p<sup>5</sup>","ORC","6.24","n/a","Greenish-yellow, disagreeable gas","1774 (Sweden)","Carl Wilhelm Scheele","Greek: chloros (greenish yellow).","17","3" },
    { "Argon","Ar","18","39.948","1.40 (@ -186Â°C)","83.8","87.3","-2","98","n/a","24.2","0.138","n/a","6.52","0.0177","85","0","1519.6","n/a","[Ne] 3s<sup>2</sup> 3p<sup>6</sup>","FCC","5.26","n/a","Colorless, tasteless, odorless noble gas","1894 (Scotland)","Sir William Ramsey, Baron Rayleigh","Greek: argos (inactive).","18","3" },
    { "Potassium","K","19","39.0983","0.856","336.8","1047","235","203","133 (+1e)","45.3","0.753","102.5","2.33","79","100","0.82","418.5","1","[Ar] 4s<sup>1</sup>","BCC","5.23","n/a","Soft, waxy, silvery-white metal","1807 (England)","Sir Humphrey Davy","English: pot ash; symbol from Latin: kalium, (alkali).","1","4" },
    { "Calcium","Ca","20","40.078","1.55","1112","1757","197","174","99 (+2e)","29.9","0.653","9.2","153.6","-201","230","1","589.4","2","[Ar] 4s<sup>2</sup>","FCC","5.58","n/a","Fairly hard, silvery-white metal","1808 (England)","Sir Humphrey Davy","Latin: calx, calcis (lime).","2","4" },
    { "Scandium","Sc","21","44.95591","2.99","1814","3104","162","144","72.3 (+3e)","15","0.556","15.8","332.7","15.8","n/a","1.36","630.8","3","[Ar] 3d<sup>1</sup> 4s<sup>2</sup>","HEX","3.31","1.594","Fairly soft, silvery-white metal","1879 (Sweden)","Lars Nilson","Latin: Scandia, Scandinavia.","3","4" },
    { "Titanium","Ti","22","47.88","4.54","1933","3560","147","132","68 (+4e) 94 (+2e)","10.6","0.523","18.8","422.6","21.9","380","1.54","657.8","4, 3","[Ar] 3d<sup>2</sup> 4s<sup>2</sup>","1.588","2.95","n/a","Shiny, dark-gray metal","1791 (England)","William Gregor","Greek: titanos (Titans).","4","4" },
    { "Vanadium","V","23","50.9415","6.11","2160","3650","134","122","59 (+5e) 74 (+3e)","8.35","0.485","17.5","460","30.7","390","1.63","650.1","5, 4, 3, 2, 0","[Ar] 3d<sup>3</sup> 4s<sup>2</sup>","BCC","3.02","n/a","Soft, ductile, silvery-white metal","1830 (Sweden)","Nils Gabriel Sefstrom","The scandinavian goddess, Vanadis.","5","4" },
    { "Chromium","Cr","24","51.9961","7.18","2130","2945","130","118","52 (+6e) 63 (+3e)","7.23","0.488","21","342","93.9","460","1.66","652.4","6, 3, 2, 0","[Ar] 3d<sup>5</sup> 4s<sup>1</sup>","BCC","2.88","n/a","Very hard, crystalline, steel-gray metal","1797 (France)","Louis Vauquelin","Greek: chroma (color).","6","4" },
    { "Manganese","Mn","25","54.93805","7.21","1517","2235","135","117","46 (+7e) 80 (+2e)","7.39","0.477","-13.4","221","-7.8","400","1.55","716.8","7, 6, 4, 3, 2, 0, -1","[Ar] 3d<sup>5</sup> 4s<sup>2</sup>","CUB","8.89","n/a","Hard, brittle, gray-white metal","1774 (Sweden)","Johann Gahn","Latin: magnes (magnet); Italian: manganese.","7","4" },
    { "Iron","Fe","26","55.847","7.874","1808","3023","126","117","64 (+3e) 74 (+2e)","7.1","0.443","13.8","~340","80.4","460","1.83","759.1","6, 3, 2, 0, -2","[Ar] 3d<sup>6</sup> 4s<sup>2</sup>","BCC","2.87","n/a","Malleable, ductile, silvery-white metal","n/a (Unknown)","Known to the ancients.","Anglo-Saxon: iron; symbol from Latin: ferrum (iron).","8","4" },
    { "Cobalt","Co","27","58.9332","8.9","1768","3143","125","116","63 (+3e) 72 (+2e)","6.7","0.456","15.48","389.1","100","385","1.88","758.1","3, 2, 0, -1","[Ar] 3d<sup>7</sup> 4s<sup>2</sup>","HEX","2.51","n/a","Hard, ductile, lustrous bluish-gray metal","1739 (Sweden)","George Brandt","German: kobold (goblin).","9","4" },
    { "Nickel","Ni","28","58.6934","8.902","1726","3005","124","115","69 (+2e)","6.6","0.443","17.61","378.6","90.9","375","1.91","736.2","3, 2, 0","[Ar] 3d<sup>8</sup> 4s<sup>2</sup>","FCC","3.52","n/a","Hard, malleable, silvery-white metal","1751 (Sweden)","Axel Cronstedt","German: kupfernickel (false copper).","10","4" },
    { "Copper","Cu","29","63.546","8.96","1356.6","2840","128","117","72 (+2e) 96 (+1e)","7.1","0.385","13.01","304.6","401","315","1.9","745","2, 1","[Ar] 3d<sup>10</sup> 4s<sup>1</sup>","FCC","3.61","n/a","Malleable, ductile, reddish-brown metal","n/a (Unknown)","Known to the ancients.","Symbol from Latin: cuprum (island of Cyprus famed for its copper mines).","11","4" },
    { "Zinc","Zn","30","65.39","7.133","692.73","1180","138","125","74 (+2e)","9.2","0.388","7.28","114.8","116","234","1.65","905.8","2","[Ar] 3d<sup>10</sup> 4s<sup>2</sup>","HEX","2.66","n/a","Bluish-silver, ductile metal","n/a (Germany)","Known to the ancients.","German: zink (German for tin).","12","4" },
    { "Gallium","Ga","31","69.723","5.91","302.93","2676","141","126","62 (+3e) 81 (+1e)","11.8","0.372","5.59","270.3","28.1","240","1.81","578.7","3","[Ar] 3d<sup>10</sup> 4s<sup>2</sup> 4p<sup>1</sup>","ORC","4.51","n/a","Soft, blue-white metal","1875 (France)","Paul-Emile Lecoq de Boisbaudran","Latin: Gallia (France).","13","4" },
    { "Germanium","Ge","32","72.61","5.323","1210.6","3103","137","122","53 (+4e) 73 (+2e)","13.6","0.322","36.8","328","60.2","360","2.01","760","4","[Ar] 3d<sup>10</sup> 4s<sup>2</sup> 4p<sup>2</sup>","DIA","5.66","n/a","Grayish-white metal","1886 (Germany)","Clemens Winkler","Latin: Germania (Germany).","14","4" },
    { "Arsenic","As","33","74.92159","5.73 (grey arsenic)","1090","876","139","120","46 (+5e) 222 (-3e)","13.1","0.328","n/a","32.4","-50.2","285","2.18","946.2","5, 3, -2","[Ar] 3d<sup>10</sup> 4s<sup>2</sup> 4p<sup>3</sup>","RHL","4.13","n/a","Steel gray, brittle semimetal","n/a (Unknown)","Known to the ancients.","Greek: arsenikon; Latin: arsenicum, (both names for yellow pigment).","15","4" },
    { "Selenium","Se","34","78.96","4.79","490","958.1","140","116","42 (+6e) 191 (-2e)","16.5","0.321 (Se-Se)","5.23","59.7","0.52","n/a","2.55","940.4","6, 4, -2","[Ar] 3d<sup>10</sup> 4s<sup>2</sup> 4p<sup>4</sup>","HEX","4.36","n/a","A soft metalloid similar to sulfur","1818 (Sweden)","Jons Jacob Berzelius","Greek: selene (moon).","16","4" },
    { "Bromine","Br","35","79.904","3.12","265.9","331.9","n/a","114","47 (+5e) 196 (-1e)","23.5","0.473 (Br-Br)","10.57 (Br-Br)","29.56 (Br-Br)","0.005","n/a","2.96","1142","7, 5, 3, 1, -1","[Ar] 3d<sup>10</sup> 4s<sup>2</sup> 4p<sup>5</sup>","ORC","6.67","n/a","Reddish-brown liquid","1826 (France)","Antoine J. Balard","Greek: bromos (stench).","17","4" },
    { "Krypton","Kr","36","83.8","2.155 (@ -153Â°C)","116.6","120.85","n/a","112","n/a","32.2","0.247","n/a","9.05","0.0095","n/a","0","1350","2","[Ar] 3d<sup>10</sup> 4s<sup>2</sup> 4p<sup>6</sup>","FCC","5.72","n/a","Dense, colorless, odorless, and tasteless gas","1898 (Great Britain)","Sir William Ramsey, M.W. Travers","Greek: kryptos (hidden).","18","4" },
    { "Rubidium","Rb","37","85.4678","1.532","312.2","961","248","216","147 (+1e)","55.9","0.36","2.2","75.8","58.2","n/a","0.82","402.8","1","[Kr] 5s<sup>1</sup>","BCC","5.59","n/a","Soft, silvery-white, highly reactive metal","1861 (Germany)","R. Bunsen, G. Kirchoff","Latin: rubidus (deep red); the color its salts impart to flames.","1","5" },
    { "Strontium","Sr","38","87.62","2.54","1042","1657","215","191","112 (+2e)","33.7","0.301","9.2","144","-35.4","n/a","0.95","549","2","[Kr] 5s<sup>2</sup>","FCC","6.08","n/a","Silvery, malleable metal","1790 (Scotland)","A. Crawford","The Scottish town, Strontian.","2","5" },
    { "Yttrium","Y","39","88.90585","4.47","1795","3611","178","162","89.3 (+3e)","19.8","0.284","11.5","367","-17.2","n/a","1.22","615.4","3","[Kr] 4d<sup>1</sup> 5s<sup>2</sup>","HEX","3.65","1.571","Silvery, ductile, fairly reactive metal","1789 (Finland)","Johann Gadolin","The Swedish village, Ytterby, where one of its minerals was first found.","3","5" },
    { "Zirconium","Zr","40","91.224","6.506","2125","4650","160","145","79 (+4e)","14.1","0.281","19.2","567","22.7","250","1.33","659.7","4","[Kr] 4d<sup>2</sup> 5s<sup>2</sup>","HEX","3.23","1.593","Gray-white, lustrous, corrosion-resistant metal","1789 (Germany)","Martin Klaproth","The mineral, zircon.","4","5" },
    { "Niobium","Nb","41","92.90638","8.57","2741","5015","146","134","69 (+5e)","10.8","0.268","26.8","680","53.7","275","1.6","663.6","5, 3","[Kr] 4d<sup>4</sup> 5s<sup>1</sup>","BCC","3.3","n/a","Shiny white, soft, ductile metal","1801 (England)","Charles Hatchet","Niobe; daughter of the mythical Greek king Tantalus.","5","5" },
    { "Molybdenum","Mo","42","95.94","10.22","2890","4885","139","130","62 (+6e) 70 (+4e)","9.4","0.251","28","~590","-138","380","2.16","684.8","6, 5, 4, 3, 2, 0","[Kr] 4d<sup>5</sup> 5s<sup>1</sup>","BCC","3.15","n/a","Silvery white, hard metal","1778 (Sweden)","Carl Wilhelm Scheele","Greek: molybdos (lead).","6","5" },
    { "Technetium","Tc","43","97.9072","11.5","2445","5150","136","127","56 (+7e)","8.5","0.243","23.8","585","50.6","n/a","1.9","702.2","7","[Kr] 4d<sup>6</sup> 5s<sup>1</sup>","HEX","2.74","1.604","Silvery-gray metal","1937 (Italy)","Carlo Perrier, Emilio Segre","Greek: technetos (artificial).","7","5" },
    { "Ruthenium","Ru","44","101.07","12.41","2583","4173","134","125","67 (+4e)","8.3","0.238","-25.5","n/a","117","n/a","2.2","710.3","8, 6, 4, 3, 2, 0, -2","[Kr] 4d<sup>7</sup> 5s<sup>1</sup>","HEX","2.7","1.584","Rare, silver-gray, extremely brittle metal","1844 (Russia)","Karl Klaus","Latin: Ruthenia (Russia).","8","5" },
    { "Rhodium","Rh","45","102.9055","12.41","2239","4000","134","125","68 (+3e)","8.3","0.244","21.8","494","150","n/a","2.28","719.5","5, 4, 3,  2, 1, 0","[Kr] 4d<sup>8</sup> 5s<sup>1</sup>","FCC","3.8","n/a","Silvery white, hard metal","1803 (England)","William Wollaston","Greek: rhodon (rose). Its salts give a rosy solution.","9","5" },
    { "Palladium","Pd","46","106.42","12.02","1825","3413","137","128","65 (+4e) 80 (+2e)","8.9","0.244","17.24","372.4","71.8","275","2.2","803.5","4,  2, 0","[Kr] 4d<sup>10</sup>","FCC","3.89","n/a","Silvery-white, soft, malleable and ductile metal","1803 (England)","William Wollaston","Named after the asteroid, Pallas, discovered in 1803.","10","5" },
    { "Silver","Ag","47","107.8682","10.5","1235.1","2485","144","134","89 (+2e) 126 (+1e)","10.3","0.237","11.95","254.1","429","215","1.93","730.5","2, 1","[Kr] 4d<sup>10</sup> 5s<sup>1</sup>","FCC","4.09","n/a","Silvery-ductile, and malleable metal","n/a (Unknown)","Known to the ancients.","Anglo-Saxon: siolful, (silver); symbol from Latin: argentium.","11","5" },
    { "Cadmium","Cd","48","112.411","8.65","594.1","1038","154","148","97 (+2e)","13.1","0.232","6.11","59.1","96.9","120","1.69","867.2","2","[Kr] 4d<sup>10</sup> 5s<sup>2</sup>","HEX","2.98","1.886","Soft, malleable, blue-white metal","1817 (Germany)","Fredrich Stromeyer","Greek: kadmeia (ancient name for calamine (zinc oxide)).","12","5" },
    { "Indium","In","49","114.818","7.31","429.32","2353","166","144","81 (+3e)","15.7","0.234","3.24","225.1","81.8","129","1.78","558","3","[Kr] 4d<sup>10</sup> 5s<sup>2</sup> 5p<sup>1</sup>","TET","4.59","n/a","Very soft, silvery-white metal","1863 (Germany)","Ferdinand Reich, T. Richter","Latin: indicum (color indigo the color it shows in a spectroscope.","13","5" },
    { "Tin","Sn","50","118.71","7.31","505.1","2543","162","141","71 (+4e) 93 (+2)","16.3","0.222","7.07","296","66.8","170","1.96","708.2","4, 2","[Kr] 4d<sup>10</sup> 5s<sup>2</sup> 5p<sup>2</sup>","TET","5.82","n/a","Silvery-white, soft, malleable and ductile metal","n/a (Unknown)","Known to the ancients.","Named after Etruscan god, Tinia; symbol from Latin: stannum (tin).","14","5" },
    { "Antimony","Sb","51","121.76","6.691","903.9","1908","159","140","62 (+6e) 245 (-3)","18.4","0.205","20.08","195.2","24.43","200","2.05","833.3","5, 3, -2","[Kr] 4d<sup>10</sup> 5s<sup>2</sup> 5p<sup>3</sup>","RHL","4.51","n/a","Hard, silvery-white, brittle semimetal","n/a (Unknown)","Known to the ancients.","Greek: anti and monos (not alone); symbol from mineral stibnite.","15","5" },
    { "Tellurium","Te","52","127.6","6.24","722.7","1263","160","136","56 (+6e) 211 (-2e)","20.5","0.201","17.91","49.8","14.3","n/a","2.1","869","6, 4, 2","[Kr] 4d<sup>10</sup> 5s<sup>2</sup> 5p<sup>4</sup>","HEX","4.45","1.33","Silvery-white, brittle semimetal","1782 (Romania)","Franz Joseph Meller von Reichenstein","Latin: tellus (earth).","16","5" },
    { "Iodine","I","53","126.90447","4.93","386.7","457.5","n/a","133","50 (+7e) 220 (-1e)","25.7","0.427 (I-I)","15.52 (I-I)","41.95 (I-I)","-0.45","n/a","2.66","1008.3","7, 5, 1, -1","[Kr] 4d<sup>10</sup> 5s<sup>2</sup> 5p<sup>5</sup>","ORC","7.72","n/a","Shiny, black nonmetallic solid","1811 (France)","Bernard Courtois","Greek: iodes (violet colored).","17","5" },
    { "Xenon","Xe","54","131.29","3.52 (@ -109Â°C)","161.3","166.1","n/a","131","n/a","42.9","0.158","n/a","12.65","0.0057","n/a","0","1170","7","[Kr] 4d<sup>10</sup> 5s<sup>2</sup> 5p<sup>6</sup>","FCC","6.2","n/a","Heavy, colorless, and odorless noble gas","1898 (England)","Sir William Ramsay; M. W. Travers","Greek: xenos (strange).","18","5" },
    { "Cesium","Cs","55","132.90543","1.873","301.6","951.6","267","235","167 (+1e)","70","0.241","2.09","68.3","35.9","n/a","0.79","375.5","1","[Xe] 6s<sup>1</sup>","BCC","6.05","n/a","Very soft, ductile, light gray metal","1860 (Germany)","Gustov Kirchoff, Robert Bunsen","Latin: coesius (sky blue); for the blue lines of its spectrum.","1","6" },
    { "Barium","Ba","56","137.327","3.5","1002","1910","222","198","134 (+2e)","39","0.192","7.66","142","-18.4","n/a","0.89","502.5","2","[Xe] 6s<sup>2</sup>","BCC","5.02","n/a","Soft, slightly malleable, silver-white metal","1808 (England)","Sir Humphrey Davy","Greek: barys (heavy or dense).","2","6" },
    { "Lanthanum","La","57","138.9055","6.15","1194","3730","187","169","101.6 (+3e)","22.5","0.197","8.5","402","13.4","132","1.1","541.1","3","[Xe] 5d<sup>1</sup> 6s<sup>2</sup>","HEX","3.75","1.619","Silvery-white, soft, malleable, and ductile metal","1839 (Sweden)","Carl Mosander","Greek: lanthanein (to be hidden).","3","9" },
    { "Cerium","Ce","58","140.115","6.757","1072","3699","181","165","92 (+4e) 103.4 (+3e)","21","0.205","5.2","398","11.3","n/a","1.12","540.1","4, 3","[Xe] 4f<sup>1</sup> 5d<sup>1</sup> 6s<sup>2</sup>","FCC","5.16","n/a","Malleable, ductile, iron-gray metal","1803 (Sweden/Germany)","W. von Hisinger, J. Berzelius, M. Klaproth","Named after the asteroid, Ceres, discovered two years before the element.","4","9" },
    { "Praseodymium","Pr","59","140.90765","6.773","1204","3785","182","165","90 (+4e) 101.3 (+3e)","20.8","0.192","11.3","331","12.5","n/a","1.13","526.6","4, 3","[Xe] 4f<sup>3</sup> 6s<sup>2</sup>","HEX","3.67","1.614","Silvery white, moderately soft, malleable, and ductile metal","1885 (Austria)","C.F. Aver von Welsbach","Greek: prasios and didymos (green twin); from its green salts.","5","9" },
    { "Neodymium","Nd","60","144.24","7.007","1294","3341","182","184","99.5 (+3e)","20.6","0.205","7.1","289","-16.5","n/a","1.14","531.5","3","[Xe] 4f<sup>4</sup> 6s<sup>2</sup>","HEX","3.66","1.614","Silvery-white, rare-earth metal that oxidizes easily in air","1925 (Austria)","C.F. Aver von Welsbach","Greek: neos and didymos (new twin).","6","9" },
//    { "Promethium","Pm","61","144.9127","7.2","1441","3000","n/a","163","97.9 (+3e)","n/a","0.185","n/a","n/a","17.9","n/a","0","536","3","[Xe] 4f<sup>5</sup> 6s<sup>2</sup>","n/a","n/a","n/a",,"1945 (United States)","J.A. Marinsky, L.E. Glendenin, C.D. Coryell","Named for the Greek god, Prometheus.","7","9" },
    { "Samarium","Sm","62","150.36","7.52","1350","2064","181","162","96.4 (+3e)","19.9","0.18","8.9","165","-13.3","166","1.17","540.1","3, 2","[Xe] 4f<sup>6</sup> 6s<sup>2</sup>","RHL","9","n/a","Silvery rare-earth metal","1853 (France)","Jean Charles Galissard de Marignac","Named after the mineral samarskite.","8","9" },
    { "Europium","Eu","63","151.965","5.243","1095","1870","199","185","95 (+3e) 109 (+2e)","28.9","0.176","n/a","176","13.9","n/a","0","546.9","3, 2","[Xe] 4f<sup>7</sup> 6s<sup>2</sup>","BCC","4.61","n/a","Soft, silvery-white metal","1901 (France)","Eugene-Antole Demarcay","Named for the continent of Europe.","9","9" },
    { "Gadolinium","Gd","64","157.25","7.9","1586","3539","179","161","93.8 (+3e)","19.9","0.23","n/a","398","-10.5","n/a","1.2","594.2","3","[Xe] 4f<sup>7</sup> 5d<sup>1</sup> 6s<sup>2</sup>","HEX","3.64","1.588","Soft, ductile, silvery-white metal","1880 (Switzerland)","Jean de Marignac","Named after the mineral gadolinite.","10","9" },
    { "Terbium","Tb","65","158.92534","8.229","1629","3296","180","159","84 (+4e) 92.3 (+3e)","19.2","0.183","n/a","389","11.1","n/a","1.2","569","4, 3","[Xe] 4f<sup>9</sup> 6s<sup>2</sup>","HEX","3.6","1.581","Soft, ductile, silvery-gray, rare-earth metal","1843 (Sweden)","Carl Mosander","Named after Ytterby, a village in Sweden.","11","9" },
    { "Dysprosium","Dy","66","162.5","8.55","1685","2835","180","159","90.8 (+3e)","19","0.173","n/a","291","10.7","n/a","n/a","567","3","[Xe] 4f<sup>10</sup> 6s<sup>2</sup>","HEX","3.59","1.573","Soft. lustrous, silvery metal","1886 (France)","Paul-Emile Lecoq de Boisbaudran","Greek: dysprositos (hard to get at).","12","9" },
    { "Holmium","Ho","67","164.93032","8.795","1747","2968","179","158","89.4 (+3e)","18.7","0.164","n/a","301","-16.2","n/a","1.23","574","3","[Xe] 4f<sup>11</sup> 6s<sup>2</sup>","HEX","3.58","1.57","Fairly soft, malleable, lustrous, silvery metal","1878 (Switzerland)","J.L. Soret","Holmia, the Latinized name for Stockholm, Sweden.","13","9" },
    { "Erbium","Er","68","167.26","9.06","1802","3136","178","157","88.1 (+3e)","18.4","0.168","n/a","317","-14.5","n/a","1.24","581","3","[Xe] 4f<sup>12</sup> 6s<sup>2</sup>","HEX","3.56","1.57","Soft, malleable, silvery metal","1843 (Sweden)","Carl Mosander","Named after the Swedish town, Ytterby.","14","9" },
    { "Thulium","Tm","69","168.93421","9.321","1818","2220","177","156","87 (+3e)","18.1","0.16","n/a","232","-16.9","n/a","1.25","589","3, 2","[Xe] 4f<sup>13</sup> 6s<sup>2</sup>","HEX","3.54","1.57","Soft, malleable, ductile, silvery metal","1879 (Sweden)","Per Theodor Cleve","Thule, ancient name of Scandinavia.","15","9" },
    { "Ytterbium","Yb","70","173.04","6.9654","1097","1466","194","n/a","85.8 (+3e) 93 (+2e)","24.8","0.145","3.35","159","-34.9","n/a","1.1","603","3, 2","[Xe] 4f<sup>14</sup> 6s<sup>2</sup>","FCC","5.49","n/a","Silvery, lustrous, malleable, and ductile metal","1878 (Switzerland)","Jean de Marignac","Named for the Swedish village of Ytterby.","16","9" },
    { "Lutetium","Lu","71","174.967","9.8404","1936","3668","175","156","85 (+3e)","17.8","0.155","n/a","414","-16.4","n/a","1.27","513","3","[Xe] 4f<sup>14</sup> 5d<sup>1</sup> 6s<sup>2</sup>","HEX","3.51","1.585","Silvery-white, hard, dense, rare-earth metal","1907 (France)","Georges Urbain","Named for the ancient name of Paris, Lutecia.","3","6" },
    { "Hafnium","Hf","72","178.49","13.31","2503","5470","167","144","78 (+4e)","13.6","0.146","-25.1","575","23","n/a","1.3","575.2","4","[Xe] 4f<sup>14</sup> 5d<sup>2</sup> 6s<sup>2</sup>","HEX","3.2","1.582","Silvery, ductile metal","1923 (Denmark)","Dirk Coster, Georg von Hevesy","Hafnia, the Latin name of Copenhagen.","4","6" },
    { "Tantalum","Ta","73","180.9479","16.654","3269","5698","149","134","68 (+5e)","10.9","0.14","24.7","758","57.5","225","1.5","760.1","5","[Xe] 4f<sup>14</sup> 5d<sup>3</sup> 6s<sup>2</sup>","BCC","3.31","n/a","Gray, heavy, hard metal","1802 (Sweden)","Anders Ekeberg","King Tantalus of Greek mythology, father of Niobe.","5","6" },
    { "Tungsten","W","74","183.84","19.3","3680","5930","141","130","62 (+6e) 70 (+4e)","9.53","0.133","-35","824","173","310","1.7","769.7","6, 5, 4, 3, 2, 0","[Xe] 4f<sup>14</sup> 5d<sup>4</sup> 6s<sup>2</sup>","BCC","3.16","n/a","Tough, steel-gray to white metal","1783 (Spain)","Juan Jose, Fausto Elhuyar","Swedish: tung sten (heavy stone): symbol from its German name wolfram.","6","6" },
    { "Rhenium","Re","75","186.207","21.02","3453","5900","137","128","53 (+7e) 72 (+4e)","8.85","0.138","34","704","48","416","1.9","759.1","5, 4, 3, 2, -1","[Xe] 4f<sup>14</sup> 5d<sup>5</sup> 6s<sup>2</sup>","HEX","2.76","1.615","Dense, silvery-white metal","1925 (Germany)","Walter Noddack, Ida Tacke, Otto Berg","Latin: Rhenus, the Rhine River.","7","6" },
    { "Osmium","Os","76","190.23","22.57","3327","5300","135","126","69 (+6e) 88 (+4e)","8.43","0.131","31.7","738","-87.6","n/a","2.2","819.8","8, 6, 4, 3, 2, 0, -2","[Xe] 4f<sup>14</sup> 5d<sup>6</sup> 6s<sup>2</sup>","HEX","2.74","1.579","Blue-white, lustrous, hard metal","1804 (England)","Smithson Tenant","Greek: osme (odor).","8","6" },
    { "Iridium","Ir","77","192.22","22.42","2683","4403","136","127","68 (+4e)","8.54","0.133","27.61","604","147","430","2.2","868.1","6, 4, 3, 2, 1, 0, -1","[Xe] 4f<sup>14</sup> 5d<sup>7</sup> 6s<sup>2</sup>","FCC","3.84","n/a","White, brittle metal","1804 (England/France)","S.Tenant, A.F.Fourcory, L.N.Vauquelin, H.V.Collet-Descoltils","Latin: iris (rainbow).","9","6" },
    { "Platinum","Pt","78","195.08","21.45","2045","4100","139","130","65 (+4e) 80 (+2e)","9.1","0.133","21.76","~470","71.6","230","2.28","868.1","4, 2, 0","[Xe] 4f<sup>14</sup> 5d<sup>9</sup> 6s<sup>1</sup>","FCC","3.92","n/a","Very heavy, soft, silvery-white metal","1735 (Italy)","Julius Scaliger","Spanish: platina (little silver).","10","6" },
    { "Gold","Au","79","196.96654","19.3","1337.58","3080","146","134","85 (+3e) 137 (+1e)","10.2","0.129","12.68","~340","318","170","2.54","889.3","3, 1","[Xe] 4f<sup>14</sup> 5d<sup>10</sup> 6s<sup>1</sup>","FCC","4.08","n/a","Soft, malleable, yellow metal","n/a (Unknown)","Known to the ancients.","Anglo-Saxon: geolo (yellow); symbol from Latin: aurum (shining dawn).","11","6" },
    { "Mercury","Hg","80","200.59","13.546 (@ +20&#176;C)","234.28","629.73","157","149","110 (+2e) 127 (+1e)","14.8","0.138","2.295","58.5","8.3","100","2","1006","2, 1","[Xe] 4f<sup>14</sup> 5d<sup>10</sup> 6s<sup>2</sup>","RHL","2.99","n/a","Heavy, silver-white metal that is in its liquid state at","n/a (Unknown)","Known to the ancients.","The Roman god Mercury; symbol from Latin: hydrargyrus (liquid silver).","12","6" },
    { "Thallium","Tl","81","204.3833","11.85","576.6","1730","171","148","95 (+3e) 147 (+1e)","17.2","0.128","4.31","162.4","46.1","96","1.62","588.9","3, 1","[Xe] 4f<sup>14</sup> 5d<sup>10</sup> 6s<sup>2</sup> 6p<sup>1</sup>","HEX","3.46","1.599","Soft, gray metal","1861 (England)","Sir William Crookes","Greek: thallos (green twig for a bright green line in its spectrum.","13","6" },
    { "Lead","Pb","82","207.2","11.35","600.65","2013","175","147","84 (+4e) 120 (+2e)","18.3","0.159","4.77","177.8","35.3","88","1.8","715.2","4, 2","[Xe] 4f<sup>14</sup> 5d<sup>10</sup> 6s<sup>2</sup> 6p<sup>2</sup>","FCC","4.95","n/a","Very soft, highly malleable and ductile, blue-white shiny metal","n/a (Unknown)","Known to the ancients.","Anglo-Saxon: lead; symbol from Latin: plumbum.","14","6" },
    { "Bismuth","Bi","83","208.98037","9.747","544.5","1883","170","146","74 (+5e) 96 (+3e)","21.3","0.124","11","172","7.9","120","2.02","702.9","5, 3","[Xe] 4f<sup>14</sup> 5d<sup>10</sup> 6s<sup>2</sup> 6p<sup>3</sup>","RHL","4.75","n/a","Hard, brittle, steel-gray metal with a pinkish tinge","n/a (Unknown)","Known to the ancients.","German: bisemutum, (white mass Now spelled wismut.","15","6" },
    { "Polonium","Po","84","208.9824","9.32","527","1235","176","146","67 (+6e)","22.7","0.125","-10","-102.9","n/a","n/a","2","813.1","6, 4, 2","[Xe] 4f<sup>14</sup> 5d<sup>10</sup> 6s<sup>2</sup> 6p<sup>4</sup>","SC","3.35","n/a","Silvery-gray metal","1898 (France/Poland)","Pierre and Marie Curie-Sklodowska","Named for Poland, native country of Marie Curie.","16","6" },
    { "Astatine","At","85","209.9871","n/a","575","610","n/a","-145","62 (+7e)","n/a","n/a","n/a","n/a","n/a","n/a","2.2","916.3","7, 5, 3, 1, -1","[Xe] 4f<sup>14</sup> 5d<sup>10</sup> 6s<sup>2</sup> 6p<sup>5</sup>","n/a","n/a","n/a","Unstable, radioactive halogen","1940 (United States)","D.R.Corson, K.R.MacKenzie, E. Segre","Greek: astatos (unstable).","17","6" },
    { "Radon","Rn","86","222.0176","4.4 (@ -62&#176;C)","202","211.4","n/a","n/a","n/a","n/a","0.094","n/a","18.1","0.0036","n/a","n/a","1036.5","n/a","[Xe] 4f<sup>14</sup> 5d<sup>10</sup> 6s<sup>2</sup> 6p<sup>6</sup>","FCC","n/a","n/a","Heavy radioactive gas","1898 (Germany)","Fredrich Ernst Dorn","Variation of the name of another element, radium.","18","6" },
    { "Francium","Fr","87","223.0197","n/a","300","950","n/a","n/a","180 (+1e)","n/a","n/a","15","n/a","n/a","n/a","0.7","~375","2","[Rn] 7s<sup>1</sup>","BCC","n/a","n/a","n/a","1939 (France)","Marguerite Derey","Named for France, the nation of its discovery.","1","7" },
    { "Radium","Ra","88","226.0254","-5.5","973","1413","n/a","n/a","143 (+2e)","45","0.12","-9.6","-113","-18.6","n/a","0.9","509","2","[Rn] 7s<sup>2</sup>","n/a","n/a","n/a","Silvery white, radioactive element","1898 (France/Poland)","Pierre and Marie Curie-Sklodowska","Latin: radius (ray).","2","7" },
    { "Actinium","Ac","89","227.0278","n/a","1320","3470","188","n/a","118 (+3e)","22.54","n/a","-10.5","-292.9","n/a","n/a","1.1","665.5","3","[Rn] 6d<sup>1</sup> 7s<sup>2</sup>","FCC","5.31","n/a","Heavy, Silvery-white metal that is very radioactive","1899 (France)","Andre-Louis Debierne","Greek: akis, aktinos (ray).","3","10" },
    { "Thorium","Th","90","232.0381","11.78","2028","5060","180","165","102 (+4e)","19.8","0.113","16.11","513.7","-54","100","1.3","670.4","4","[Rn] 6d<sup>2</sup> 7s<sup>2</sup>","FCC","5.08","n/a","Gray, soft, malleable, ductile, radioactive metal","1828 (Sweden)","Jons Jacob Berzelius","Named for Thor, Norse god of thunder.","4","10" },
    { "Protactinium","Pa","91","231.03588","15.37","2113","4300","161","n/a","89 (+5e) 113 (+3e)","15","0.121","16.7","481.2","n/a","n/a","1.5","n/a","5, 4","[Rn] 5f<sup>2</sup> 6d<sup>1</sup> 7s<sup>2</sup>","TET","3.92","n/a","Silvery-white, radioactive metal","1917 (England/France)","Fredrich Soddy, John Cranston, Otto Hahn, Lise Meitner","Greek: proto and actinium (parent of actinium); it forms actinium when it radioactively decays.","5","10" },
    { "Uranium","U","92","238.0289","19.05","1405.5","4018","138","142","80 (+6e) 97 (+4e)","12.5","0.115","12.6","417","27.5","n/a","1.38","686.4","6, 5, 4, 3","[Rn] 5f<sup>3</sup> 6d<sup>1</sup> 7s<sup>2</sup>","ORC","2.85","n/a","Silvery-white, dense, ductile and malleable, radioactive metal.","1789 (Germany)","Martin Klaproth","Named for the planet Uranus.","6","10" },
    { "Neptunium","Np","93","237.048","20.25","913","4175","130","n/a","95 (+4e) 110 (+3e)","21.1","n/a","-9.6","336","-6.3","n/a","1.36","n/a","6, 5, 4, 3","[Rn] 5f<sup>4</sup> 6d<sup>1</sup> 7s<sup>2</sup>","ORC","4.72","n/a","Silvery metal","1940 (United States)","E.M. McMillan, P.H. Abelson","Named for the planet Neptune.","7","10" },
    { "Plutonium","Pu","94","244.0642","19.84","914","3505","151","n/a","93 (+4e) 108 (+3e)","n/a","n/a","2.8","343.5","-6.7","n/a","1.28","491.9","6, 5, 4, 3","[Rn] 5f<sup>6</sup> 7s<sup>2</sup>","MCL","n/a","n/a","Silvery-white, radioactive metal","1940 (United States)","G.T.Seaborg, J.W.Kennedy, E.M.McMillan, A.C.Wohl","Named for the planet Pluto.","8","10" },
    { "Americium","Am","95","243.0614","13.67","1267","2880","173","n/a","92 (+4e) 107 (+3e)","20.8","n/a","-10","238.5","n/a","n/a","1.3","n/a","6, 5, 4, 3","[Rn] 5f<sup>7</sup> 7s<sup>2</sup>","n/a","n/a","n/a","Silvery-white, radioactive metal","1945 (United States)","G.T.Seaborg, R.A.James, L.O.Morgan, A.Ghiorso","Named for the American continent, by analogy with europium.","9","10" },
    { "Curium","Cm","96","247.0703","13.51","1340","n/a","299","n/a","n/a","18.28","n/a","n/a","n/a","n/a","n/a","1.3","-580","4, 3","[Rn] 5f<sup>7</sup> 6d<sup>1</sup> 7s<sup>2</sup>","n/a","n/a","n/a","Silvery, malleable, synthetic radioactive metal","1944 (United States)","G.T.Seaborg, R.A.James, A.Ghiorso","Named in honor of Pierre and Marie Curie.","10","10" },
    { "Berkelium","Bk","97","247.0703","13.25","n/a","n/a","297","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","1.3","-600","4, 3","[Rn] 5f<sup>9</sup> 7s<sup>2</sup>","n/a","n/a","n/a","Radionactive synthetic metal","1949 (United States)","G.T.Seaborg, S.G.Tompson, A.Ghiorso","Named after Berkeley, California the city of its discovery.","11","10" },
    { "Californium","Cf","98","251.0796","15.1","900","n/a","295","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","1.3","-610","4, 3","[Rn] 5f<sup>10</sup> 7s<sup>2</sup>","n/a","n/a","n/a","Powerful neutron emitter","1950 (United States)","G.T.Seaborg, S.G.Tompson, A.Ghiorso, K.Street Jr.","Named after the state and University of California.","12","10" },
    { "Einsteinium","Es","99","252.083","n/a","n/a","1130","292","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","1.3","-620","3","[Rn] 5f<sup>11</sup> 7s<sup>2</sup>","n/a","n/a","n/a","Radioactive, synthetic metal","1952 (United States)","Argonne, Los Alamos, U of Calif","Named in honor of the scientist Albert Einstein.","13","10" },
    { "Fermium","Fm","100","257.0951","n/a","1800","n/a","290","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","1.3","-630","3","[Rn] 5f<sup>12</sup> 7s<sup>2</sup>","n/a","n/a","n/a","Radioactive, synthetic metal","1953 (United States)","Argonne, Los Alamos, U of Calif","Named in honor of the scientist Enrico Fermi.","14","10" },
    { "Mendelevium","Md","101","258.1","n/a","1100","n/a","287","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","1.3","-635","3","[Rn] 5f<sup>13</sup> 7s<sup>2</sup>","n/a","n/a","n/a","Radioactive, synthetic metal","1955 (United States)","G.T.Seaborg, S.G.Tompson, A.Ghiorso, K.Street Jr.","Named in honor of the scientist Dmitri Ivanovitch Mendeleyev, who devised the periodic table.","15","10" },
    { "Nobelium","No","102","259.1009","n/a","1100","n/a","285","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","1.3","-640","3, 2","[Rn] 5f<sup>14</sup> 7s<sup>2</sup>","n/a","n/a","n/a","Radioactive, synthetic metal.","1957 (Sweden)","Nobel Institute for Physics","Named in honor of Alfred Nobel, who invented dynamite and founded Nobel prize.","16","10" },
    { "Lawrencium","Lr","103","262.11","n/a","n/a","n/a","282","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","3","[Rn] 5f<sup>14</sup> 6d<sup>1</sup> 7s<sup>2</sup>","n/a","n/a","n/a","Radioactive, synthetic metal","1961 (United States)","A.Ghiorso, T.Sikkeland, A.E.Larsh, R.M.Latimer","Named in honor of Ernest O. Lawrence, inventor of the cyclotron.","3","7" },
    { "Rutherfordium","Rf","104","[261]","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","[Rn] 5f<sup>14</sup> 6d<sup>2</sup> 7s<sup>2</sup>","n/a","n/a","n/a","Radioactive, synthetic metal","1969 (United States)","A. Ghiorso, et al","Named in honor of Ernest Rutherford","4","7" },
    { "Dubnium","Db","105","[262]","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","[Rn] 5f<sup>14</sup> 6d<sup>3</sup> 7s<sup>2</sup>","n/a","n/a","n/a","Radioactive, synthetic metal","1970 (United States)","A. Ghiorso, et al","Named in honor of Otto Hahn","5","7" },
    { "Seaborgium","Sg","106","[266]","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","[Rn] 5f<sup>14</sup> 6d<sup>4</sup> 7s<sup>2</sup>","n/a","n/a","n/a","Radioactive, synthetic metal","1974 (USSR/United States)","Soviet Nuclear Research/ U. of Cal at Berkeley","Named in honor of Glenn Seaborg, American physical chemist known for research on transuranium elements.","6","7" },
    { "Bohrium","Bh","107","[264]","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","[Rn] 5f<sup>14</sup> 6d<sup>5</sup> 7s<sup>2</sup>","n/a","n/a","n/a","n/a","Radioactive, synthetic metal","1976 (Germany)","Heavy Ion Research Laboratory (HIRL)","Named in honor of Niels Bohr","7","7" },
    { "Hassium","Hs","108","[269]","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","[Rn] 5f<sup>14</sup> 6d<sup>6</sup> 7s<sup>2</sup>","n/a","n/a","n/a","n/a","1984 (Germany)","Heavy Ion Research Laboratory (HIRL)","Named in honor of Henri Hess, Swiss born Russian chemist known for work in thermodydamics.","8","7" },
    { "Meitnerium","Mt","109","[268]","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","[Rn] 5f<sup>14</sup> 6d<sup>7</sup> 7s<sup>2</sup>","n/a","n/a","n/a","n/a","1982 (Germany)","Heavy Ion Research Laboratory (HIRL)","Named in honor of Lise Mietner","9","7" },
    { "Ununnilium","Uun","110","[269]","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","[Rn] 5f<sup>14</sup> 6d<sup>9</sup> 7s<sup>1</sup>","n/a","n/a","n/a","n/a","1994 (Germany)","Heavy Ion Research Laboratory (HIRL)","Un (one) nun (one) nilium (zero)","10","7" },
    { "Unununium","Uuu","111","[272]","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","[Rn] 5f<sup>14</sup> 6d<sup>10</sup> 7s<sup>1</sup>","n/a","n/a","n/a","n/a","1994 (Germany)","Heavy Ion Research Laboratory (HIRL)","Un (one) nun (one) unium (one)","11","7" },
    { "Ununbium","Uub","112","[277]","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","[Rn] 5f<sup>14</sup> 6d<sup>10</sup> 7s<sup>2</sup>","n/a","n/a","n/a","n/a","1996 (n/a)","n/a","n/a","12","7" },
    { "Ununtrium","Uut","113","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","[Rn] 5f<sup>14</sup> 6d<sup>10</sup> 7s<sup>2</sup> 7p<sup>1</sup>","n/a","n/a","n/a","n/a","n/a (n/a)","n/a","n/a","13","7" },
    { "Ununquadium","Uuq","114","[289]","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","[Rn] 5f<sup>14</sup> 6d<sup>10</sup> 7s<sup>2</sup> 7p<sup>2</sup>","n/a","n/a","n/a","n/a","1999 (n/a)","n/a","n/a","14","7" },
    { "Ununpentium","Uup","115","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","[Rn] 5f<sup>14</sup> 6d<sup>10</sup> 7s<sup>2</sup> 7p<sup>3</sup>","n/a","n/a","n/a","n/a","n/a (n/a)","n/a","n/a","15","7" },
    { "Ununhexium","Uuh","116","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","[Rn] 5f<sup>14</sup> 6d<sup>10</sup> 7s<sup>2</sup> 7p<sup>4</sup>","n/a","n/a","n/a","n/a","1999 (n/a)","n/a","n/a","16","7" },
    { "Ununseptium","Uus","117","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","[Rn] 5f<sup>14</sup> 6d<sup>10</sup> 7s<sup>2</sup> 7p<sup>5</sup>","n/a","n/a","n/a","n/a","n/a (n/a)","n/a","n/a","17","7" },
    { "Ununoctium","Uuo","118","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","[Rn] 5f<sup>14</sup> 6d<sup>10</sup> 7s<sup>2</sup> 7p<sup>6</sup>","n/a","n/a","n/a","n/a","1999 (n/a)","n/a","n/a","18","7" }
    };

    static std::string cleanStringForNum(std::string str)
    {
        std::string cleanStr = str;
        cleanStr.replace(cleanStr.find("("), cleanStr.find(")")-cleanStr.find("("), "");
        return cleanStr;
    }

    static double getAtomicWeight(std::string element_symbol)
    {
        int colIDX = 3;
        if( strncmp((headers[colIDX]).c_str(), "Atomic Weight", strlen("Atomic Weight")) )
        {
            ConsoleOutput::WriteWarning("Probably accessing wrong values\n"
                    "ColIDX : " +std::to_string(colIDX)+ " Header Entry : " + (headers[colIDX]),
                    thisclassname, "getAtomicWeight()");
        }
        for (unsigned int row = 0; row < periodicTable.size(); row++)
        {
            if(strncmp(periodicTable[row][1].c_str(), element_symbol.c_str(), strlen(element_symbol.c_str())) )
            {
                return atof(cleanStringForNum(periodicTable[row][colIDX]).c_str());
                break;
            }
        }
    }

    static double getAtomicVolume(std::string element_symbol)
    {
        int colIDX = 10;
        if( strncmp((headers[colIDX]).c_str(), "Atomic Volume", strlen("Atomic Volume")) )
        {
            ConsoleOutput::WriteWarning("Probably accessing wrong values\n"
                    "ColIDX : " +std::to_string(colIDX)+ " Header Entry : " + (headers[colIDX]),
                    thisclassname, "getAtomicVolume()");
        }
        for (unsigned int row = 0; row < periodicTable.size(); row++)
        {
            if(strncmp(periodicTable[row][1].c_str(), element_symbol.c_str(), strlen(element_symbol.c_str())) )
            {
                return atof(cleanStringForNum(periodicTable[row][colIDX]).c_str());
                break;
            }
        }
    }

    static double getDensity(std::string element_symbol)
    {
        int colIDX = 4;
        if( strncmp((headers[colIDX]).c_str(), "Density", strlen("Density")) )
        {
            ConsoleOutput::WriteWarning("Probably accessing wrong values\n"
                    "ColIDX : " +std::to_string(colIDX)+ " Header Entry : " + (headers[colIDX]),
                    thisclassname, "getAtomicVolume()");
        }
        for (unsigned int row = 0; row < periodicTable.size(); row++)
        {
            if(strncmp(periodicTable[row][1].c_str(), element_symbol.c_str(), strlen(element_symbol.c_str())) )
            {
                return atof(cleanStringForNum(periodicTable[row][colIDX]).c_str());
                break;
            }
        }
    }

    //TODO convert to struct
 //    struct Element
 //    {
 //        std::string "name";
 //        std::string "symbol";
 //        std::string "atomic_number";
 //        std::string "atomic_weight";
 //        std::string "density";
 //        std::string "melting_point";
 //        std::string "boiling_point";
 //        std::string "atomic_radius";
 //        std::string "covalent_radius";
 //        std::string "ionic_radius";
 //        std::string "atomic_volume";
 //        std::string "specific_heat";
 //        std::string "fusion_heat";
 //        std::string "evaporation_heat";
 //        std::string "termal_conductivity";
 //        std::string "debye_temperature";
 //        std::string "pauling_negativity_number";
 //        std::string "first_ionizing_energy";
 //        std::string "oxidation_states";
 //        std::string "electronic_configuration";
 //        std::string "lattice_structure";
 //        std::string "lattice_constant";
 //        std::string "lattice_c-a_ratio";
 //        std::string "info_appearance";
 //        std::string "info_discovery_date";
 //        std::string "info_discovered_by";
 //        std::string "info_named_after";
 //        std::string "pos_x";
 //        std::string "pos_y"
 //    };

protected:
private:
};

}  // namespace openphase

#endif /* IELEMENTS_H */
