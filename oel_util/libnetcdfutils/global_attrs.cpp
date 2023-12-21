#include "global_attrs.h"
#include <ctime>
#include <netcdf>
#include <libgen.h>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

string call_sequence(int argc, char* argv[]) {

  // get processing timestamp, to the second
  time_t tnow = time(nullptr);
  char prodtime[80];
  strftime(prodtime, 80, "%Y-%m-%dT%XZ", gmtime(&tnow));

  // append calling sequence
  string callseq = prodtime;
  callseq.append(":");
  for (int i=0; i<argc; i++) {
    callseq.append(" ");
    callseq.append(basename(argv[i]));
  }
  return callseq;
}

string get_history(NcFile *ncfile) {
  string history;
  NcGroupAtt att = ncfile->getAtt("history");
  if (!att.isNull()) att.getValues(history);
  if (history.length() > 0) { history.append("; "); }
  return history;
}

void set_global_attrs(string filename,
                      string history,
                      string doi) {

  // Open file for append
  NcFile *outfile = new NcFile(filename, NcFile::write);

  // set standard global metadata
  outfile->putAtt("institution", "NASA Goddard Space Flight Center, Ocean Biology Processing Group");
  outfile->putAtt("license", "https://science.nasa.gov/earth-science/earth-science-data/data-information-policy/");
  outfile->putAtt("project", "Ocean Biology Processing Group");
  outfile->putAtt("naming_authority", "gov.nasa.gsfc.oceandata");

  outfile->putAtt("creator_name", "NASA/GSFC/OBPG");
  outfile->putAtt("creator_email", "data@oceancolor.gsfc.nasa.gov");
  outfile->putAtt("creator_url", "https://oceancolor.gsfc.nasa.gov");
  outfile->putAtt("publisher_name", "NASA/GSFC/OB.DAAC");
  outfile->putAtt("publisher_email", "data@oceancolor.gsfc.nasa.gov");
  outfile->putAtt("publisher_url", "https://oceancolor.gsfc.nasa.gov");

  outfile->putAtt("keywords_vocabulary", "NASA Global Change Master Directory (GCMD) Science Keywords");
  outfile->putAtt("Conventions", "CF-1.8 ACDD-1.3");
  outfile->putAtt("standard_name_vocabulary", "CF Standard Name Table v79");

  // set doi info
  if (!doi.empty()) {
    outfile->putAtt("identifier_product_doi_authority", "https://dx.doi.org");
    outfile->putAtt("identifier_product_doi", doi);
  }

  // set processing history
  if (!history.empty()) {
    outfile->putAtt("history", history);
  }

  // set date_created
  outfile->putAtt("date_created", unix2isodate(now(), 'G'));

  // done
  outfile->close();
}
