#ifndef PERIODICTABLE_H
#define PERIODICTABLE_H

const int kZmax=118;

char element_symbol[kZmax][4]={
{"H"},{"He"},{"Li"},{"Be"},{"B"},{"C"},{"N"},{"O"},
{"F"},{"Ne"},{"Na"},{"Mg"},{"Al"},{"Si"},{"P"},{"S"},
{"Cl"},{"Ar"},{"K"},{"Ca"},{"Sc"},{"Ti"},{"V"},{"Cr"},
{"Mn"},{"Fe"},{"Co"},{"Ni"},{"Cu"},{"Zn"},{"Ga"},{"Ge"},
{"As"},{"Se"},{"Br"},{"Kr"},{"Rb"},{"Sr"},{"Y"},{"Zr"},
{"Nb"},{"Mo"},{"Tc"},{"Ru"},{"Rh"},{"Pd"},{"Ag"},{"Cd"},
{"In"},{"Sn"},{"Sb"},{"Te"},{"I"},{"Xe"},{"Cs"},{"Ba"},
{"La"},{"Ce"},{"Pr"},{"Nd"},{"Pm"},{"Sm"},{"Eu"},{"Gd"},
{"Tb"},{"Dy"},{"Ho"},{"Er"},{"Tm"},{"Yb"},{"Lu"},{"Hf"},
{"Ta"},{"W"},{"Re"},{"Os"},{"Ir"},{"Pt"},{"Au"},{"Hg"},
{"Tl"},{"Pb"},{"Bi"},{"Po"},{"At"},{"Rn"},{"Fr"},{"Ra"},
{"Ac"},{"Th"},{"Pa"},{"U"},{"Np"},{"Pu"},{"Am"},{"Cm"},
{"Bk"},{"Cf"},{"Es"},{"Fm"},{"Md"},{"No"},{"Lr"},{"Rf"},
{"Db"},{"Sg"},{"Bh"},{"Hs"},{"Mt"},{"Ds"},{"Rg"},{"Uub"},
{"Uut"},{"Uuq"},{"Uup"},{"Uuh"},{"Uus"},{"Uuo"}
};

enum element_enum {elem_H,elem_He,elem_Li,elem_Be,elem_B,elem_C,elem_N,elem_O,
		   elem_F,elem_Ne,elem_Na,elem_Mg,elem_Al,elem_Si,elem_P,elem_S,
		   elem_Cl,elem_Ar,elem_K,elem_Ca,elem_Sc,elem_Ti,elem_V,elem_Cr,
		   elem_Mn,elem_Fe,elem_Co,elem_Ni,elem_Cu,elem_Zn,elem_Ga,elem_Ge,
		   elem_As,elem_Se,elem_Br,elem_Kr,elem_Rb,elem_Sr,elem_Y,elem_Zr,
		   elem_Nb,elem_Mo,elem_Tc,elem_Ru,elem_Rh,elem_Pd,elem_Ag,elem_Cd,
		   elem_In,elem_Sn,elem_Sb,elem_Te,elem_I,elem_Xe,elem_Cs,elem_Ba,
		   elem_La,elem_Ce,elem_Pr,elem_Nd,elem_Pm,elem_Sm,elem_Eu,elem_Gd,
		   elem_Tb,elem_Dy,elem_Ho,elem_Er,elem_Tm,elem_Yb,elem_Lu,elem_Hf,
		   elem_Ta,elem_W,elem_Re,elem_Os,elem_Ir,elem_Pt,elem_Au,elem_Hg,
		   elem_Tl,elem_Pb,elem_Bi,elem_Po,elem_At,elem_Rn,elem_Fr,elem_Ra,
		   elem_Ac,elem_Th,elem_Pa,elem_U,elem_Np,elem_Pu,elem_Am,elem_Cm,
		   elem_Bk,elem_Cf,elem_Es,elem_Fm,elem_Md,elem_No,elem_Lr,elem_Rf,
		   elem_Db,elem_Sg,elem_Bh,elem_Hs,elem_Mt,elem_Ds,elem_Rg,elem_Uub,
		   elem_Uut,elem_Uuq,elem_Uup,elem_Uuh,elem_Uus,elem_Uuo};

// Most common isotope 
// Incomplete: Make sure you complete for elements you 
// want to use.
int a_stable[kZmax]={
  1,0,0,0,11,12,14,16,
  19,0,23,24,27,28,0,0,
  35,40,0,0,0,48,0,0,
  0,56,59,0,64,65,0,0,
  0,0,0,0,0,0,0,0,
  93,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,207,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0
};

char element_name[kZmax][15]=
{
{"Hydrogen"},{"Helium"},{"Lithium"},{"Beryllium"},{"Boron"},{"Carbon"},
{"Nitrogen"},{"Oxygen"},{"Fluorine"},{"Neon"},{"Sodium"},{"Magnesium"},
{"Aluminium"},{"Silicon"},{"Phosphorus"},{"Sulphur"},{"Chlorine"},{"Argon"},
{"Potassium"},{"Calcium"},{"Scandium"},{"Titanium"},{"Vanadium"},{"Chromium"},
{"Manganese"},{"Iron"},{"Cobalt"},{"Nickel"},{"Copper"},{"Zinc"},{"Gallium"},
{"Germanium"},{"Arsenic"},{"Selenium"},{"Bromine"},{"Krypton"},{"Rubidium"},
{"Strontium"},{"Yttrium"},{"Zirconium"},{"Niobium"},{"Molybdenum"},{"Techenium"},
{"Ruthenium"},{"Rhodium"},{"Palladium"},{"Silver"},{"Cadmium"},{"Indium"},{"Tin"},
{"Antimony"},{"Tellurium"},{"Iodine"},{"Xenon"},{"Cesium"},{"Barium"},{"Lanthanum"},
{"Cerium"},{"Praseodymium"},{"Neodymium"},{"Prometium"},{"Samarium"},{"Europium"},
{"Gadolinium"},{"Terbium"},{"Dysprosium"},{"Holmium"},{"Erbium"},{"Thulium"},{"Ytterbium"}
,{"Lutetium"},{"Hafnium"},{"Tantalum"},{"Tungsten"},{"Rhenium"},{"Osmium"},{"Iridium"},
{"Platinum"},{"Gold"},{"Mercury"},{"Thallium"},{"Lead"},{"Bismuth"},{"Polonium"},{"Astatine"}
,{"Radon"},{"Francium"},{"Radium"},{"Actinium"},{"Thorium"},{"Protactinium"},{"Uranium"},
{"Neptunium"},{"Plutonium"},{"Americium"},{"Curium"},{"Berkelium"},{"Californium"},{"Einsteinium"},
{"Fermium"},{"Mendelevium"},{"Nobelium"},{"Lawrencium"},{"Rutherfordium"},{"Dubinium"},
{"Seaborgium"},{"Bohrium"},{"Hassium"},{"Meitnerium"},{"Darmstadtium"},{"Roentgenium"},
{"Ununbium"},{"Ununtrium"},{"Ununquadium"},{"Ununpentium"},{"Ununhexium"},{"Ununseptium"},
{"Ununoctium"}
};

#endif
