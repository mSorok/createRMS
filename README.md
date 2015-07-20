## CREATION OF REACTION MOLECULAR SIGNATURES FROM METACYC REACTIONS
===================================================================
Last modification : 20/07/2015


__Requirements:__

* MetaCyc flat files (to download)
  * reactions.dat
  * MOL files in the "MetaCyc-MOLfiles" directory
* ChemAxon's molconvert
* molsig (http://molsig.sourceforge.net/)
* a MySQL database

----------------------------------------------------------------------------------------------------------------------

### Step 1: Protonation et aromatization of MOL files


`$ python add_hydrogens_aromatization.py molconvert_PATH`


----------------------------------------------------------------------------------------------------------------------

### Step 2: Compute Molecular Signatures on all molecules

`$ python molsigLauncher.py MolFiles_FULL_aroma sscan n`

* arg1: directory where are located the MOLfiles for which compute the molecular signature
* arg2: *scan* / *sscan* / *fsscan* molsig parameter
* arg3: *o*/*n* = with/without aromatization


----------------------------------------------------------------------------------------------------------------------

### Step 3: Compute the Reaction Molecular Signatures


`$ python rms_compute.py sscan n > rms_sscan.txt`

* arg1: *scan* / *sscan* / *fsscan* molsig parameter
* arg2: *o*/*n* = with/without aromatization

----------------------------------------------------------------------------------------------------------------------

### Step 4: Insertion of RMS in database

    
#### DB tables creaction

```MySQL
    
 -- Reaction - RMS correpondency table
DROP TABLE IF EXISTS Reaction_RMS_CPD;
CREATE TABLE Reaction_RMS_CPD(
height INT(3),
MR_id VARCHAR(255),
RMS TEXT(80000)
);
ALTER TABLE Reaction_RMS_CPD ADD INDEX(MR_id);


    
-- RMS - MD5 encoded RMS (RMSh)
DROP TABLE IF EXISTS RMS_RMSh_CPD;
CREATE TABLE RMS_RMSh_CPD(
RMS TEXT(80000),
RMSh VARCHAR(255)
);
ALTER TABLE RMS_RMSh_CPD ADD INDEX (RMSh);
    

-- Reaction - RMSh correspondency table
DROP TABLE IF EXISTS Reaction_RMSh_CPD;
CREATE TABLE Reaction_RMSh_CPD(
MR_id VARCHAR(255),
RMSh VARCHAR(255),
height INT(3),
reaction_type enum('balanced','unbalanced','no-pwy') DEFAULT 'unbalanced'
);
ALTER IGNORE TABLE Reaction_RMSh_CPD ADD UNIQUE INDEX(MR_id, RMSh, height);
    
```	
	
	

     
#### Insert file in database
* insert `rms_sscan.txt` in `Reaction_RMS_CPD` table 
        

#### Transform text RMS in MD5 encoded RMSh

```MySQL
INSERT IGNORE INTO RMS_RMSh_CPD
SELECT RMS, MD5(RMS)
FROM Reaction_RMS_CPD
GROUP BY RMS;
     
INSERT IGNORE INTO Reaction_RMSf_CPD(MR_id, RMSh,height)
SELECT MR_id, RMSh, height
FROM Reaction_RMS_CPD 
    INNER JOIN RMS_RMSh_CPD USING(RMS);
        

    
-- AS at height 0 balanced reactions have their RMS = '0.0', it's a way to detect unbalanced reactions
UPDATE Reaction_RMSh_CPD AS RRF, Reaction_RMS_CPD AS RR
SET RRF.reaction_type = 'balanced'
WHERE RRF.MR_id = RR.MR_id AND RR.height=0 AND RR.RMS = '0.0'; 
    
     
    
DELETE FROM Reaction_RMSf_CPD WHERE RMSh IS NULL;
```

----------------------------------------------------------------------------------------------------------------------

###Step 5: Creation of intelligible RMS identifiers: RMSid

_Extract from the database in `MR_RMSf_chain.txt` using this SQL request:_

```MySQL
SELECT MR_id, GROUP_CONCAT(RMSf ORDER BY diameter SEPARATOR '$') 
FROM Reaction_RMSf_CPD WHERE reaction_type = 'balanced' 
GROUP BY MR_id ORDER BY MR_id DESC ;
```


_Launch of the RMSid creator_
`$ python rms_id_creator.py MR_RMSf_chain.txt > rmsf_d_rmsid_cpd.txt`



_Integration of RMSids in the database_

```MySQL
DROP TABLE IF EXISTS RMSh_RMSid_CPD;
CREATE TABLE RMSh_RMSid_CPD(
 RMSh VARCHAR(255),
 height INT(11),
 RMSid VARCHAR(255)
);
ALTER TABLE RMSh_RMSid_CPD ADD UNIQUE INDEX (RMSh,RMSid);
```


* insert `rmsf_d_rmsid_cpd.txt` in the table `RMSf_RMSid_CPD`


```MySQL
DROP TABLE IF EXISTS Reaction_RMSid_CPD;
CREATE TABLE Reaction_RMSid_CPD(
 MR_id VARCHAR(255),
 RMSid VARCHAR(255),
 height INT(3)
);
INSERT INTO Reaction_RMSid_CPD
SELECT t1.MR_id, t2.RMSid,t2.height
FROM Reaction_RMSh_CPD AS t1
  INNER JOIN RMSh_RMSid_CPD AS t2 USING (RMSh,height);
```

----------------------------------------------------------------------------------------------------------------------


