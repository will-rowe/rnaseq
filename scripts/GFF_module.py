#!/usr/bin/env python
__author__ = 'Will Rowe'
__mail__ = "will.rowe@liverpool.ac.uk"

import os
import sys
import time
import gzip # used to handle gzipped files
import collections # used for OrderedDict
import getpass # for identifying username
import logging # for appending to the pipeline log file
logger = logging.getLogger(__name__)


####
# Information
####
""" The GFF module is used in the Hinton Lab RNA-seq pipeline

This module is required to:

    * check annotation files for compliance with GFF3 formatting
    * create new annotation files (in GFF3 format)
    * perform TPM calculations
    * generate reports for single samples

Background:

    User supplied annotation files (often in Excel worksheets) frequently break the Hinton Lab RNA-seq pipeline and checks are needed to prevent this.

Extra information:

GFF (generic feature format) is a standard file format for storing genomic features in a text file
GFF files are plain text, 9 column, tab-delimited files
See -> https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
"""


####
# Classes
####
class GFF_annotationFile:
    """ Class summary: creates an object from a GFF annotation file.

    notes on handling GFF attributes:
    * currently all possible GFF attributes are included in the fileFeautres class attribute (inside a nested dictionary) - if they are not present in the file, they are added with 'no_value' set
    * when parsing the attributes from column 9, inverted commas are not removed (see GFF3 documentation)
    * trailing semi-colon is removed from final entry if it is present
    """
    # set the GFF_annotationFile fileCount class variable (value shared among all instances of this class - accessed as GFF_annotationFile.fileCount)
    file_count = 0


    def __init__(self, gff_filename):
        """
        set the class constructor with input parameters
        """
        GFF_annotationFile.file_count += 1
        self.gff_filename = os.path.abspath(gff_filename)
        self.file_feature_count = 0
        self.file_errors_dict = {}
        self.file_info = []


    def reportError(self):
        """
        simple handler for minor errors - currently results in sys.exit but could be modified to issue warnings etc.
        """
        logging.error('The GFF module encountered an error')
        logging.error('errors found:\t{}' .format(len(self.file_errors_dict)))
        for key in self.file_errors_dict:
            logging.error('error:\t{}\tcause:\t{}' .format(self.file_errors_dict[key], key))
        sys.exit(1)


    def run_gffRead(self):
        """
        read in the GFF file, perform basic checks and populate GFF_annotationFile object attributes
        """
        self.file_features_dict = collections.OrderedDict()
        self.file_fasta_dict = {}

        # check file exists
        if not os.path.isfile(self.gff_filename):
            self.file_errors_dict[self.gff_filename] = 'supplied GFF file does not exist --> create file before reading . . .'
            self.reportError()

        # allow for gz compression
        openFunc = gzip.open if self.gff_filename.endswith(".gz") else open

        # read in the file, assess each line and run appropriate methods to populate object attributes
        with openFunc(self.gff_filename) as gff:
            for line in gff:
                line = line.rstrip('\n')

                # ignore empty lines
                if line == '':
                    continue

                # read GFF file comments
                if line.startswith('#'):
                    self.file_info.append('## comment line from original GFF file: {}' .format(line))

                # read feature entries, making sure lines have correct number of columns (9)
                elif not line.startswith('#') and len(line.split('\t')) == 9:
                    self.run_featureUnpack(line)
                    self.file_feature_count += 1

                # read fasta entires (doesn't account for N's etc. at start of fasta sequence)
                elif line.startswith('>'):
                    fastaHeader = line[1:]
                elif line.upper().startswith(('A','T','C','G')) and fastaHeader:
                    if fastaHeader in self.file_fasta_dict.keys():
                        self.file_fasta_dict[fastaHeader].append(line)
                    else:
                        self.file_fasta_dict[fastaHeader] = [line]

                # create error if line hasn't already been handled
                else:
                    self.file_errors_dict[line] = 'line not handled'

        # once all lines have been read, check if no features found in file
        if self.file_feature_count == 0:
            self.file_errors_dict[line] = 'no features identified in file'

        # if any errors were found whilst reading file --> report and exit
        if self.file_errors_dict:
            self.reportError()


    def run_featureUnpack(self, line):
        """
        unpacks each gff feature into a list (plus an ordered dictionary of attributes). The list is then saved as a dictionary entry for the gff object.
        """
        seqid, source, seqtype, start, end, score, strand, phase, attributes_list = line.split('\t')
        attributes_dictionary = collections.OrderedDict()

        # get values for feature attributes. the line below: 1. removes trailing ';' from attributes_list, splits the string on every ';' occurence, splits each of these substrings into key and value around first "=" encountered
        try:
            for attribute in attributes_list.rstrip(';').split(';'):
                 key, value = attribute.split('=',1)
                 attributes_dictionary[key] = value
        except:
            self.file_errors_dict[line] = 'bad attribute formatting'

        # use the ID value (should be unique - we'll check at the end of the unpack) as the key for the whole feature in the self.file_features_dict dictionary
        if "ID" in attributes_dictionary.keys():
            attribute_ID = attributes_dictionary["ID"]
        elif not "ID" in attributes_dictionary.keys():
            attributes_dictionary["ID"] = 'no_value'
        if not "Name" in attributes_dictionary.keys():
            attributes_dictionary["Name"] = 'no_value'
        if not "Alias" in attributes_dictionary.keys():
            attributes_dictionary["Alias"] = 'no_value'
        if not "Parent" in attributes_dictionary.keys():
            attributes_dictionary["Parent"] = 'no_value'
        if not "Target" in attributes_dictionary.keys():
            attributes_dictionary["Target"] = 'no_value'
        if not "Gap" in attributes_dictionary.keys():
            attributes_dictionary["Gap"] = 'no_value'
        if not "Derives_from" in attributes_dictionary.keys():
            attributes_dictionary["Derives_from"] = 'no_value'
        if not "Note" in attributes_dictionary.keys():
            attributes_dictionary["Note"] = 'no_value'
        if not "Dbxref" in attributes_dictionary.keys():
            attributes_dictionary["Dbxref"] = 'no_value'
        if not "Ontology_term" in attributes_dictionary.keys():
            attributes_dictionary["Ontology_term"] = 'no_value'
        if not "Is_circular" in attributes_dictionary.keys():
            attributes_dictionary["Is_circular"] = 'no_value'
        feature_list = seqid, source, seqtype, start, end, score, strand, phase, attributes_dictionary

        # check the feature doesn't already exist in the GFF file and then add the feature to the feature dictionary
        if attribute_ID in self.file_features_dict.keys():
            self.file_errors_dict[line] = 'duplicate attribute ID found'
        else:
            self.file_features_dict[attribute_ID] = feature_list


    def run_gffCheck(self):
        """
        run sanity checks on an annotation file
        at the moment this kills the pipeline, these checks could be handled as warnings instead
        """
        # make sure the GFF file has been read
        if self.file_feature_count == 0:
            self.run_gffRead()

        # count the number of features (use an arbitary check for now)
        if self.file_feature_count <= 1000:
            self.file_errors_dict[self.file_feature_count] = 'fewer than 1000 features found in the GFF file (this is an arbitary threshold) - only:'

        # are CDS and ncRNA present (use the seqtype for each feature)?
        if sum(self.file_features_dict[feature][2].upper() == "CDS" for feature in self.file_features_dict) == 0 or sum(self.file_features_dict[feature][2].upper() == "NCRNA" for feature in self.file_features_dict) == 0:
            self.file_errors_dict['NOcdsORrna'] = 'no CDS or ncRNA features found in the GFF file'

        # check fasta entries are present
        if not self.file_fasta_dict:
            self.file_errors_dict['NOfasta'] = 'no fasta entries found in the GFF file'

        # if any errors were found whilst running checks --> report and exit
        if self.file_errors_dict:
            self.reportError()

        self.file_info.append('## file checks: original GFF file ({}) passed pipeline checks for formatting and content' .format(self.gff_filename))


    def make_gffFile(self, output_filename, write_fasta=True):
        """
        create an output GFF file from a GFF_annotationFile object
        the original feature order should have been maintained thanks to the orderedDict BUT if we have added new features they will be at the end of the feature list regardless of start value
        """
        # make sure the GFF has the appropriate attributes filled and passes the GFF checks
        if not self.file_feature_count:
            self.run_gffCheck()

        # format the object so it's ready for writing
        feature_list = []
        for key in self.file_features_dict:
            # process the first 8 fields (single items), more processing required for 9th field (dictionary of attributes)
            feature_line = '{}\t' .format('\t'.join(map(str, self.file_features_dict[key][0:8])))

            # format attributes dictionary and add to the end of the feature_line
            for attribute_name in self.file_features_dict[key][8]:
                feature_line += '{}={};' .format(attribute_name, self.file_features_dict[key][8][attribute_name])

            # add the formatted feature line to the list of features
            feature_list.append(feature_line)

        # get the user who is writing the new annotation file
        user = getpass.getuser()

        with open(output_filename, "w") as gff_output_file:
            gff_output_file.write('## GFF version 3\n## file created by: {}\n## file created on: {}\n' .format(user, time.strftime('%d %b %Y @ %H:%M:%S')))
            if self.file_info:
                gff_output_file.write('{}\n' .format('\n'.join(map(str, self.file_info))))
            gff_output_file.write('{}' .format('\n'.join(map(str, feature_list))))

            # write fasta to end of file (this will be unordered)
            if write_fasta == True:
                for key in self.file_fasta_dict:
                    gff_output_file.write('\n>{}\n{}\n' .format(key, ''.join(map(str, self.file_fasta_dict[key]))))


    def get_featureLength(self, attribute_query):
        """
        method to calculate length of each feature ---> for use in the TPM calculations etc.
        """
        length = int(self.file_features_dict[attribute_query][4]) - int(self.file_features_dict[attribute_query][3])
        return length

    def get_fastaSeqs(self):
        fasta_seqs = ''
        for key in self.file_fasta_dict:
            fasta_seqs += '>{}\n{}\n' .format(key, ''.join(map(str, self.file_fasta_dict[key])))
        return fasta_seqs




class sampleReport(GFF_annotationFile):
    """ Class summary: a subclass of GFF_annotationFile that is used to hold sample data
    The class needs to have the reference GFF data that was used to align sample reads etc.


    This class builds on the standard GFF annotation file by adding count and TPM values to each feature
    """


    def __init__(self, gff_filename, count_data_file):
        """
        set the class constructor with input parameters
        """
        GFF_annotationFile.__init__(self, gff_filename)
        self.count_data_file = count_data_file
        self.count_values_dict = {}
        self.tpm_values_dict = {}

        # make sure the sampleReport has had the reference GFF loaded in
        if not self.file_feature_count:
            self.run_gffCheck()


    def get_featureCounts(self, fileformat='featureCounts'):
        """
        adds a count value to the count_data_file dictionary for each feature
        can support count files generated by either featureCounts or htseq-count (if specified)
        """
        # check file exists
        if not os.path.isfile(self.count_data_file):
            self.file_errors_dict[self.count_data_file] = 'supplied count file does not exist . . .'
            self.reportError()

        # allow for gz compression
        openFunc = gzip.open if self.gff_filename.endswith(".gz") else open

    	# open the count file and IF feature in our annotation, add each count to the self.count_values_dict using the feature name as the key
        with openFunc(self.count_data_file) as counts:
            if fileformat == 'featureCounts':
                for line in counts:
                    # ignore comments/headers
                    if line.startswith('#'):
                        continue
                    if line.startswith('Geneid'):
                        continue
                    entry = line.rstrip('\n').split('\t')
                    if not len(entry) == 7:
                        continue
                    try:
                        feature, count = entry[0], int(entry[6])
                    except:
                        self.file_errors_dict[line] = 'can\'t process count data file - is it featureCounts output?'
                        self.reportError()
                    if feature not in self.count_values_dict.keys() and feature in self.file_features_dict.keys():
                        self.count_values_dict[feature] = count
                    elif feature not in self.file_features_dict.keys():
                        self.file_errors_dict[line] = 'feature in count file is NOT present in the reference GFF . . .'
                        self.reportError()
                    else:
                        self.file_errors_dict[line] = 'can\'t process count data file - possible duplicate or formatting issue'
                        self.reportError()
            elif fileformat == 'htseq':
                for line in counts:
                    # don't include extra information from htseq (eg. __no_feature) in the count dictionary
                    if feature.startswith("__"):
                        continue
                    line = line.split("\t")
                    try:
                        feature, count = line[0], int(line[1])
                    except:
                        self.file_errors_dict[line] = 'can\'t process count data file - is it htseq-count output?'
                        self.reportError()
                else:
                    self.file_errors_dict[line] = 'can\'t process count data file - format not supported'
                    self.reportError()

                if feature not in self.count_values_dict.keys() and feature in self.file_features_dict.keys():
                    self.count_values_dict[feature] = count
                elif feature not in self.file_features_dict.keys():
                    self.file_errors_dict[line] = 'feature in count file is NOT present in the reference GFF . . .'
                    self.reportError()
                else:
                    self.file_errors_dict[line] = 'can\'t process count data file - possible duplicate or formatting issue'
                    self.reportError()


    def get_tpmValues(self):
        """
        adds a TPM value to the tpm_values_dict dictionary for each feature
        """
        # make sure the count data has been read
        if not self.count_values_dict:
            self.get_featureCounts()

        # calculate the sum of all rates for the sample
        sum_of_rates = 0
        for feature in self.count_values_dict.keys():
            sum_of_rates += float(self.count_values_dict[feature]) / float(self.get_featureLength(feature))

        # calculate TPM value per feature
        for feature in self.count_values_dict.keys():
            tpm_value = (float(self.count_values_dict[feature]) / float(self.get_featureLength(feature))) * (1.0 / sum_of_rates) * (10**6)
            self.tpm_values_dict[feature] = tpm_value


    def make_sampleReport(self, output_filename):
        """
        create a report for a single sample
        contains the feature ID, a description, the location, count values and TPM values
        """
        # make sure TPMs have been calculated
        if not self.tpm_values_dict:
            self.get_tpmValues()

        # format the object so it's ready for writing
        feature_list = []
        for key in self.count_values_dict:
            feature_line = '{}\t{}\t{}\t{}\t{}' .format(key, self.file_features_dict[key][3], self.file_features_dict[key][4], self.count_values_dict[key], self.tpm_values_dict[key])
            feature_list.append(feature_line)

        # get the user who is writing the report
        user = getpass.getuser()

        with open(output_filename, "w") as sample_report:
            sample_report.write('## TPM report\n## file created by: {}\n## file created on: {}\n' .format(user, time.strftime('%d %b %Y @ %H:%M:%S')))
            if self.file_info:
                sample_report.write('{}\n' .format('\n'.join(map(str, self.file_info))))
            sample_report.write('feature\tstart\tend\tcount value\ttpm value\n')
            sample_report.write('{}' .format('\n'.join(map(str, feature_list))))
