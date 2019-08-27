# Cross-project collaboration

Let's take an example of two different projects, Project A: a group of researchers at Radboud working with complex genetics and Project B: a group of researchers at VUMC working with clinical genetics. Both projects need to download the same data from the UK Biobank resource and perform their own analysis. Here is a high level overview of a proposed solution on Spider:

![Spider scientific catalog](/extras/Spider_scientific_catalog.png)

## Roles 

Each project in this example has its own environment where the project members can collaborate by sharing data, software or workflows. The project team members are organized with different roles and mandates:

* Data manager(s): designated data dissemination manager; responsible for the management of project-owned data
* Software manager(s): designated software manager; responsible to install and maintain the project-owned software
* Normal users: researchers who focus on their data analysis and visualization
* External collaborators: researchers who don't have access to the platform but they can be given direct access (through a web link) in selected directories containing designated data products

Each project may have also some project-specific data on external storage managed by a different group that can be accessed through standard protocols and APIs (e.g. webdav, S3, sftp). 

Let's say now that the normal user groups of both project A Complex genetics and project B Clinical genetics need access to the same data. Instead of downloading multiple copies of the same data on the same platform, we keep a single copy UK Biobank data and arrange read-only access for both groups if they comply with the data regulations. Here, the `catalog manager` is responsible for populating the catalog and deciding which project groups have read access to the catalog.

## Permissions

The table below summarizes the permissions on the different project directories:

| Directories vs. Access Roles | /project/ProjectA/Data | /project/ProjectA/Software | /project/ProjectA/Share | /project/ProjectA/Public | /project/home/ProjectA-user | /catalog/UK-biobank |
| -------------------------------|---|---|---|---|---|---|
| Project A Data manager         |rwx|r-x|rwx|rwx|---|---|    
| Project A Software manager     |r-x|rwx|rwx|rwx|---|---| 
| Project A normal user          |r-x|r-x|rwx|rwx|rwx|r--|  
| Project B normal user          |---|---|---|r--|---|r--|  
| Catalog manager                |---|---|---|r--|---|rwx| 



> **_Food for brain:_**
>
> * Can you find any existing catalogs on the platform that your project `spidercourse` has access to? 
> * You can share data with other collaborators through different project folders: `project/share`, `project/public`, `/catalog`. Can you think of scanarios that best suit each option?





