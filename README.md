######################################################################################
#      CREATION OF REACTION MOLECULAR SIGNATURES FROM METACYC REACTIONS              #
######################################################################################

Last modification : 11/06/2015


Requirements:
	- MetaCyc flat files
		- reactions.dat
		- MetaCyc-MOLfiles
	- ChemAxon's molconvert
	- molsig (http://molsig.sourceforge.net/)
	- a MySQL database

***********************************************************************************************************************



    
***********************************************************************************************************************

Protonation et aromatization of MOL files


$ python add_hydrogens_aromatization.py molconvert_PATH


***********************************************************************************************************************
Launch molsig on all MOLfiles

	1st arg : MOLfiles directory to use
	2nd arg : molsig type (scan, sscan, fsscan)
	3nd arg : o/n with/without aromatization



$ python molsigLauncher.py MolFiles_FULL_aroma sscan n


***********************************************************************************************************************
Compute RMS

	1st arg : molsig type (scan, sscan, fsscan)
	2nd arg : o/n with/without aromatization
	
$ python rms_compute.py sscan n > rms_sscan.txt

***********************************************************************************************************************

Creation and insertion of RMS in database

    
TABLES CREATION
    
    -- Reaction - RMS correpondency table
    DROP TABLE IF EXISTS Reaction_RMS_CPD;
    CREATE TABLE Reaction_RMS_CPD(
    diameter INT(3),
    MR_id VARCHAR(255),
    RMS TEXT(80000)
    );
    
    ALTER TABLE Reaction_RMS_CPD ADD INDEX(MR_id);


    
    -- RMS - MD5 encoded RMS (RMSf)
    DROP TABLE IF EXISTS RMS_RMSf_CPD;
    CREATE TABLE RMS_RMSf_CPD(
    RMS TEXT(80000),
    RMSf VARCHAR(255)
    );
    ALTER TABLE RMS_RMSf_CPD ADD INDEX (RMSf);
    

    -- Reaction - RMSf correspondency table
    DROP TABLE IF EXISTS Reaction_RMSf_CPD;
    CREATE TABLE Reaction_RMSf_CPD(
    MR_id VARCHAR(255),
    RMSf VARCHAR(255),
    diameter INT(3),
    reaction_type enum('balanced','unbalanced','no-pwy') DEFAULT 'unbalanced'
    );
    ALTER IGNORE TABLE Reaction_RMSf_CPD ADD UNIQUE INDEX(MR_id, RMSf,diameter);
	
	
	

     
INSERTION IN DABABASE
	o insert rms_sscan.txt in Reaction_RMS_CPD table 
        

Transform text RMS in MD5 encoded RMSf


    INSERT IGNORE INTO RMS_RMSf_CPD
    SELECT RMS, MD5(RMS)
    FROM Reaction_RMS_CPD
    GROUP BY RMS;
     
    INSERT IGNORE INTO Reaction_RMSf_CPD(MR_id, RMSf,diameter)
    SELECT MR_id, RMSf, diameter
    FROM Reaction_RMS_CPD 
        INNER JOIN RMS_RMSf_CPD USING(RMS);
        

    
    -- AS at height 0 balanced reactions have their RMS = '0.0', it's a way to detect unbalanced
	reactions
    UPDATE Reaction_RMSf_CPD AS RRF, Reaction_RMS_CPD AS RR
    SET RRF.reaction_type = 'balanced'
    WHERE RRF.MR_id = RR.MR_id AND RR.diameter=0 AND RR.RMS = '0.0'; 
    
     
    
    DELETE FROM Reaction_RMSf_CPD WHERE RMSf IS NULL;

***********************************************************************************************************************
RMSid creation 

Extract from the database in MR_RMSf_chain.txt using this SQL request:

	SELECT MR_id, GROUP_CONCAT(RMSf ORDER BY diameter SEPARATOR '$') 
	FROM Reaction_RMSf_CPD WHERE reaction_type = 'balanced' 
	GROUP BY MR_id ORDER BY MR_id DESC ;




$ python rms_id_creator.py MR_RMSf_chain.txt > rmsf_d_rmsid_cpd.txt





Integration of RMSids in the database


	DROP TABLE IF EXISTS RMSf_RMSid_CPD;
	CREATE TABLE RMSf_RMSid_CPD(
    	RMSf VARCHAR(255),
    	diameter INT(11),
    	RMSid VARCHAR(255)
	);
	ALTER TABLE RMSf_RMSid_CPD ADD UNIQUE INDEX (RMSf,RMSid);



Insert rmsf_d_rmsid_cpd.txt in table RMSf_RMSid_CPD





	DROP TABLE IF EXISTS Reaction_RMSid_CPD;
	CREATE TABLE Reaction_RMSid_CPD(
    	MR_id VARCHAR(255),
    	RMSid VARCHAR(255),
    	diameter INT(3)
	);
	INSERT INTO Reaction_RMSid_CPD
	SELECT t1.MR_id, t2.RMSid,t2.diameter
	FROM Reaction_RMSf_CPD AS t1
    	INNER JOIN RMSf_RMSid_CPD AS t2 USING (RMSf,diameter);



***********************************************************************************************************************


