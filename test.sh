#
# Remove the currently installed eco and install the version for test
#
R CMD remove eco

#
# Install HJ08003/eco for the test
#
echo "===================== test HJ08003/eco  ============================"
Rscript -e 'devtools::install_github("HJ08003/eco", force = TRUE)'
rm -f hj.txt; Rscript test.R > hj.txt


#                                                                                                                                                      
# Remove the current eco
#
R CMD remove eco

#
# Install kosukeimai/eco for the test
#
echo "===================== test kosukeimai/eco =========================="
Rscript -e 'devtools::install_github("kosukeimai/eco", force = TRUE)'
rm -f ki.txt; Rscript test.R > ki.txt

#
# Display the difference
#
echo ""
echo ""
echo "====== The difference between Github eco versions of KI and HJ ======"
echo ""
diff ki.txt hj.txt
echo "===================== Above is the difference  ======================"

#
# Make sure always install the official version at the end of the test
#
# R CMD remove eco
# Rscript -e 'devtools::install_github("kosukeimai/eco", force = TRUE)'
# 
