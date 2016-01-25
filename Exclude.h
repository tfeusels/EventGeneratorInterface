#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoShape.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>

using std::string;



void RecursiveExclude(string gOptRootGeomTopVol, TGeoVolume *volume, bool exclude)
{
  if (!volume) return; // just a precaution
  
  // check if the "current volume" is in the gOptRootGeomTopVol list
  if (volume->GetName()) // non-null pointer to volume name?
    {
      const char *name = volume->GetName();
      size_t length = strlen(name);
      if (length) // non-empty volume name?
        {
	  size_t ind = 0;
	  while ( (ind = gOptRootGeomTopVol.find_first_of("+-", ind)) != std::string::npos )
            {
	      ind += 1;
	      if (ind == gOptRootGeomTopVol.length()) break; // just a precaution
	      if ( (!(gOptRootGeomTopVol.compare(ind, length, name))) &&
		   ( ((ind + length) == gOptRootGeomTopVol.length()) ||
		     (gOptRootGeomTopVol[(ind + length)] == '+') ||
		     (gOptRootGeomTopVol[(ind + length)] == '-') ||
		     (gOptRootGeomTopVol[(ind + length)] == ' ') ) )
                {
		  // a "match" is found ... now check what to do with it
		  if ( gOptRootGeomTopVol[(ind - 1)] == '+')
		    exclude = false;
		  else if ( gOptRootGeomTopVol[(ind - 1)] == '-')
                        exclude = true;
		  break; // we are done
                }
            }
        }
    }
  
#if defined(DEBUG_RECURSIVE_EXCLUDE)
  if(!exclude){
    std::cout << volume->GetName()
              << " <" << volume->GetMedium()->GetName() << ">"
              << " : " << exclude << " :";
  }
#endif /* defined(DEBUG_RECURSIVE_EXCLUDE) */
    
    // "exclude" the "current volume" if requested
    if (exclude)
      {
        static TGeoMaterial *matVacuum = ((TGeoMaterial *)0);
        static TGeoMedium   *Vacuum = ((TGeoMedium *)0);
        
        if (!Vacuum)
	  {
#if defined(DEBUG_RECURSIVE_EXCLUDE)
            std::cout << " Creating the Vaccum material and medium :";
#endif /* defined(DEBUG_RECURSIVE_EXCLUDE) */
            // Actually ... one should check if the "Vacuum" TGeoMaterial and
            // TGeoMedium are already defined in the geometry and, if found,
            // re-use them ... but I was too lazy to implement it here, sorry.
            if (!matVacuum) matVacuum = new TGeoMaterial("Vacuum", 0.0, 0.0, 0.0);
            if (matVacuum) Vacuum = new TGeoMedium("Vacuum", 1, matVacuum);
	  }
        
        if (Vacuum) volume->SetMedium(Vacuum); // "exclude" volume
      }
    
#if defined(DEBUG_RECURSIVE_EXCLUDE)
    std::cout << " <" << volume->GetMedium()->GetName() << ">"
              << std::endl;
#endif /* defined(DEBUG_RECURSIVE_EXCLUDE) */
    
    // proceed with all daughters of the "current volume"
    Int_t nd = volume->GetNdaughters();
    for (Int_t i = 0; i < nd; i++)
      {
        if (volume->GetNode(i)) // non-null pointer to node?
	  RecursiveExclude(gOptRootGeomTopVol, volume->GetNode(i)->GetVolume(), exclude);
      }
}

void PrepareVolume(string gOptRootGeomTopVol, TGeoVolume* topvol){
  std::cout << gOptRootGeomTopVol << std::endl;
   if ( (gOptRootGeomTopVol[0] == '+') || (gOptRootGeomTopVol[0] == '-') ) {

#if defined(DEBUG_RECURSIVE_EXCLUDE)
     //        std::cout << "set in main ... rgeom->SetTopVolName(\"\"); // \"master volume\""
     //            << std::endl;
#endif /* defined(DEBUG_RECURSIVE_EXCLUDE) */
        // "exclude" (make them out of "Vacuum") all non-selected volumes,
        // the first "action" character in gOptRootGeomTopVol decides ...
        RecursiveExclude(gOptRootGeomTopVol, topvol, (gOptRootGeomTopVol[0] == '+'));
   } else {
     std::cout << "gOptRootGeomTopVol does not begin with +/- : "
	       << gOptRootGeomTopVol << std::endl;
#if defined(DEBUG_RECURSIVE_EXCLUDE)
     std::cout << "set in main ... rgeom->SetTopVolName(\""
	       << gOptRootGeomTopVol << "\"); // user selected volume"
	       << std::endl;
#endif /* defined(DEBUG_RECURSIVE_EXCLUDE) */
   }
   
   return; 
}



//____________________________________________________________________________

/* From file testRecursiveExclude.cxx by Jacek M. Holeczek */
