import shutil
import fileinput
import json
import itertools



class HeaderRenaming:

    FILL_LEN = 7
    HEADER_PREFIX = 'C'

    def rename_fasta_headers(self, fasta_filepath, fasta_filepath_new=None, 
                             cache_filepath=None):

        """
        if `fasta_filepath_new` not specified, will inplace edit `fasta_filepath`
        Will save mapping in file `self.cache_filepath` if is given, else will save mapping in self.new2old
        """
        new2old = {}
        # work on `fasta_filepath_new`
        if fasta_filepath_new:
            shutil.copyfile(fasta_filepath, fasta_filepath_new)
        else:
            fasta_filepath_new = fasta_filepath
        #
        for i, line in zip(itertools.count(), fileinput.input(fasta_filepath_new, inplace=True)):
            if line.startswith('>'):
                header_old = line.strip()[1:]
                header_new = self.HEADER_PREFIX + str(i).zfill(self.FILL_LEN)
                new2old[header_new] = header_old
                print('>' + header_new)
            else:
                print(line, end='')
        # mapping saving mode
        if cache_filepath:
            with open(cache_filepath, 'w') as fp:
                json.dump(new2old, fp)
            self.cache_filepath = cache_filepath
            self.new2old = None
        else:
            self.cache_mapping_filepth = None
            self.new2old = new2old



    def revert_fasta_headers(self, fasta_filepath, fasta_filepath_new=None):
        # work on `fasta_filepath_new`
        if fasta_filepath_new:
            shutil.copyfile(fasta_filepath_old, fasta_filepath_new)
        else:
            fasta_filepath_new = fasta_filepath
        # mapping saving mode
        if self.cache_filepath:
            with open(self.cache_filepath) as fp:
                new2old = json.load(fp)
        else:
            new2old = self.new2old
        #
        for line in fileinput.input(fasta_filepath_new, inplace=True):
            if line.startswith('>'):
                header_new = line.strip()[1:]
                header_old = new2old[header_new]
                print('>' + header_old)
            else:
                print(line, end='')
        







