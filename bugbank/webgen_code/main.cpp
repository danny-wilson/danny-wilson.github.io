/*	webgen.cpp
 *	This program reads from webgen.txt two columns, one a list of proto-html
 *	filenames and one a list of titles for those pages. A template is then
 *	applied to generate complete html files.
 *  This program also reads navigation.txt to build the navigation bar
 *  and footer.txt to build the footer.
 *  A style sheet style.css is assumed to exist in the same folder.
 *
 *	Required associated files:
 *		webgen.txt
 *		navigation.txt
 *      footer.txt
 *      style.css (at download)
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>

using namespace std;

void error(const char* errTxt);

int main(int argc, char* argv[]) {
	ifstream webgen("webgen.txt");
	if(!webgen || !webgen.is_open()) error("webgen.txt not found");
	
	vector<string> filename;
	vector<string> title;

	const char* tab = "\t";
	string temp;
	while(!webgen.eof()) {
		getline(webgen,temp);
		int gap = temp.find_first_of(tab);
		if(gap!=temp.find_last_of(tab)) {
			string errTxt = "Found multiple tabs in line " + temp;
			error(errTxt.c_str());
		}
		if(gap==-1) {
			if(webgen.eof()) break;
			string errTxt = "No tab found in line " + temp;
			error(errTxt.c_str());
		}
		filename.push_back(temp.substr(0,gap));
		title.push_back(temp.substr(gap+1,-1));
	}
	int i,j,n;
	if((n=filename.size())!=title.size()) error("Number of file names doesn't match number of titles");

	for(i=0;i<n;i++) {
		int dot = filename[i].find_last_of(".");
		string htmlname = filename[i].substr(0,dot) + ".html";
		if(htmlname==filename[i]) {
			string errTxt = "Source and target filename the same " + filename[i];
			error(errTxt.c_str());
		}
		/* Open prototype navigation bar */
		ifstream navig("navigation.txt");
		if(!navig || !navig.is_open()) error("navigation.txt not found");
		/* Open prototype content */
		ifstream raw1(filename[i].c_str());
		if(!raw1 || !raw1.is_open()) error(filename[i].c_str());
		/* Open prototype footer */
		ifstream footer("footer.txt");
		if(!footer || !footer.is_open()) error("footer.txt not found");

		/* No longer actually a frame, this program builds a stand-alone html page*/
		ofstream frame(htmlname.c_str());
		if(!frame || !frame.is_open()) error(htmlname.c_str());
		/* HTML header */
		frame << "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\" \"http://www.w3.org/TR/html4/loose.dtd\">" << endl;
		frame << "<html>" << endl;
		frame << "<head>" << endl;
		frame << tab << "<title>" << title[i] << "</title>" << endl;
		frame << tab << "<meta http-equiv=\"Content-Type\" content=\"text/html; charset=iso-8859-1\">" << endl;
		/* Load the style sheet */
		frame << tab << "<link rel=\"StyleSheet\" href=\"style.css\" type=\"text/css\" media=\"screen\">" << endl;
		/* JavaScript: if framed break out of frame */
		frame << tab << "<SCRIPT LANGUAGE=\"JAVASCRIPT\" TYPE=\"TEXT/JAVASCRIPT\">" << endl;
		frame << tab << "<!--" << endl;
		frame << tab << "if(window != top) top.location.href = location.href" << endl;
		frame << tab << "//-->" << endl;
		frame << tab << "</SCRIPT>" << endl;
		frame << "</head>" << endl;
		/* Body */
		frame << "<body>" << endl;
		/* Title banner */
		frame << tab << "<div id=\"title\">" << endl;
		frame << tab << tab << "<div id=\"header\">" << endl;
		frame << tab << tab << tab << "<h1>Bugbank</h1>" << endl;
		frame << tab << tab << "</div>" << endl;
		frame << tab << "</div>" << endl;
		/* Navigation bar */
		frame << tab << "<div id=\"navigation\" class=\"class1\">" << endl;
		while(!navig.eof()) {
			/* Copy in the contents of the navigation bar */
			getline(navig,temp);
			frame << tab << tab << temp << endl;
		}
		frame << tab << "</div>" << endl;
		/* Main content */
		frame << tab << "<div id=\"content\">" << endl;
		while(!raw1.eof()) {
			/* Copy in the contents of the main page */
			getline(raw1,temp);
			frame << tab << tab << temp << endl;
		}
		frame << tab << "</div>" << endl;
		/* Footer */
		frame << tab << "<div id=\"footer\">" << endl;
		while(!footer.eof()) {
			/* Copy in the contents of the footer */
			getline(footer,temp);
			frame << tab << tab << temp << endl;
		}
		frame << tab << "</div>" << endl;
		/* Stat-counter code */
		frame << tab << "<!-- Start of StatCounter Code -->" << endl;
		frame << tab << "<script type=\"text/javascript\" language=\"javascript\">" << endl;
		frame << tab << "<!-- " << endl;
		frame << tab << "var sc_project=1989169; " << endl;
		frame << tab << "var sc_invisible=1; " << endl;
		frame << tab << "var sc_partition=18; " << endl;
		frame << tab << "var sc_security=\"e3f4bbdd\"; " << endl;
		frame << tab << "//-->" << endl;
		frame << tab << "</script>" << endl;
		frame << tab << endl;
		frame << tab << "<script type=\"text/javascript\" language=\"javascript\" src=\"http://www.statcounter.com/counter/frames.js\"></script><noscript><a href=\"http://www.statcounter.com/\" target=\"_blank\"><img  src=\"http://c19.statcounter.com/counter.php?sc_project=1989169&java=0&security=e3f4bbdd&invisible=1\" alt=\"counter\" border=\"0\"></a> </noscript>" << endl;
		frame << tab << "<!-- End of StatCounter Code -->" << endl;
		/***/
		/* Google Analytics code */
		frame << tab << "<script type=\"text/javascript\">" << endl;
		frame << tab << "var gaJsHost = ((\"https:\" == document.location.protocol) ? \"https://ssl.\" : \"http://www.\");" << endl;
		frame << tab << "document.write(unescape(\"%3Cscript src='\" + gaJsHost + \"google-analytics.com/ga.js'" << endl;
		frame << tab << "type='text/javascript'%3E%3C/script%3E\"));" << endl;
		frame << tab << "</script>" << endl;
		frame << tab << "<script type=\"text/javascript\">" << endl;
		frame << tab << "var pageTracker = _gat._getTracker(\"UA-5003457-1\");" << endl;
		frame << tab << "pageTracker._initData();" << endl;
		frame << tab << "pageTracker._trackPageview();" << endl;
		frame << tab << "</script>" << endl;
		/***/
		frame << "</body>" << endl;
		frame << "</html>" << endl;
		frame.close();
		footer.close();
		raw1.close();
		navig.close();
	}

	webgen.close();
	return 0;
}

void error(const char* errTxt) {
	cout << "ERROR: " << errTxt << endl;
	exit(13);
}
