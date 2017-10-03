# maize-fractionation

See other [README](doc/README.md) for full details.

# Fractionation Project

This project mostly follows the guidelines from the R package -- ProjectTemplate

Here is the architecture: 

- projects/ 
	- fractionation/
		- cache/
		- data/
		- diagnostics/
		- doc/
		- graphs/
		- munge/
		- reports/
		- src/  
		README.md  
		TODO.md  
		.gitignore  

Here is summary of what these directories contain:

**cache**: Here you’ll store any data sets that (a) are generated during a preprocessing step and (b) don’t need to be regenerated every single time you analyze your data.  

**data**: Here you'll store your raw data files.  

**diagnostics**: Here you can store any scripts you use to diagnose your data sets for corruption or problematic data points.  

**doc**: Here you can store any documentation that you've written about your analysis.

**graphs**: Here you can store any graphs that you produce.  

**munge**: Here you can store any preprocessing or data munging code for your project. For example, if you need to add columns at runtime, merge normalized data sets or globally censor any data points, that code should be store in the munge directory. 

**reports**: Here you can store any output reports, such as HTML or LaTex version of tables, that you produce.

**src**: Here you'll store your final, robust analysis scripts. 

**README**: In this file, you should write some notes to help orient any newcomers to your project.

**TODO**: In this file, you should write a list of future improvements and bug fixes that you plant ot make to your analyses.
