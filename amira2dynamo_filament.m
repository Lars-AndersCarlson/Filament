function table_out = amira2dynamo_filament(points,point_IDs,offsetx,offsety,offsetz,apix,ypix,tomo_nr,offset_reg_nr,template_fil,model_path,cropdist,crop_dphi)
%function table_out = amira2dynamo_filament(points,point_IDs,offsetx,offsety,offsetz,apix,ypix,tomo_nr,offset_reg_nr,template_fil,model_path,cropdist,crop_dphi)
%
%This function imports data from Amira's filament tracing to Dynamo and
%saves the data in the form of model files. It also returns a table
%containing all the crop points.
%
%The data from Amira should first be saved as a spreadsheet. Then the X, Y
%and Z coordinates in the tab "Points" should be saved as one csv file, and
%the column "Point IDs" from the tab Segments should be saved as another
%csv file (in LibreOffice, the text delimiter must be changed to '). 
%
%NOTE that the point ID csv file needs to be further edited to remove all
%commas before the first number on each line and also remove the single
%quote (') before the first and after the last number. Otherwise (for some
%mysterious reason) the last number in each filament (i.e. line of the
%matrix) is set to -1.
%
%Then import each csv into MATLAB. For each, the right area has to be 
%marked in the import wizard. For the point IDs, the "unimportable values" 
%must be replaced by -1.
%
%INPUT:
%points         The matrix that contains the coordinates of each point, as
%               imported from the csv.
%point_IDs      The matrix that contains the point IDs, as imported from the
%               csv.
%offsetx        These are the values displayed as "from" when clicking on the
%               volume used for filament tracing in Amira. (first value)
%offsety        (second value)
%offsetz        (third value)
%apix           pixel size in Angstrom of the unbinned tomogram that will be
%               used by Dynamo for cropping.
%ypix           Size in Y (second dimension) of the unbinned tomogram as per
%               above (3710 in the case of an unbinned K2 tomogram).
%tomo_nr        The tomogram number that will be entered into column 20 of the
%               resulting table.
%offset_reg_nr  Each filament gets a unique number in column 21 in the
%               resulting table. The lowest number will be offset_reg_nr + 1
%template_fil   Path to a model file that will be used as a template for
%               the new model files. This should be a filamentWithTorsion
%               model and should contain the right info in the cvolume
%               structure related to the volume that all the filaments come
%               from. Easiest is just to manually save one filament in dcm.
%model_path     The directory in which the filament models should be saved.
%cropdist       The distance in unbinned pixels at which subvolumes should
%               be cropped.
%crop_dphi      The angle dphi used for the cropping, in degrees
%
%OUTPUT:
%table_out      The combined table with all the crop points for all
%               filaments.
%
%usage example:
% table_out=amira2dynamo_filament(points,pointIDs,274.624,-274.624,4393.98,17.164/4,3710,10,0,'filaments_1/tomograms/volume_10/models/template/template.omd','filaments_1/tomograms/volume_10/models/',10,60);
% dwrite(table_out,'filaments_1/tomograms/volume_10/croptable.tbl');
% dtcrop lars_filaments_100px_2.Data/indices_column20.doc filaments_1/tomograms/volume_10/croptable.tbl vol10_amiratest.Data 100
% daverage('vol10_amiratest.Data','table','filaments_1/tomograms/volume_10/croptable.tbl','fcompensate',1,'o','filaments_1/tomograms/volume_10/vol10_testaverage.em');
% 
%2020-11-09 LAC

nfil=size(point_IDs,1); %number of filaments
disp(['There will be ' num2str(nfil) ' filaments']) %debug output

template_filament=dread(template_fil);

offsetx=offsetx/apix;%in pixels
offsety=offsety/apix;%in pixels
offsetz=offsetz/apix;%in pixels

newpoints=zeros(size(points,1),3);%these transformations go from Amira to 3dmod/dynamo coordinate system (i.e. to pixels)
newpoints(:,1)=points(:,1)/apix-offsetx;
newpoints(:,2)=ypix-(points(:,2)/apix-offsety);
newpoints(:,3)=points(:,3)/apix-offsetz;

%size(newpoints) %debug output
%newpoints %debug output
%disp(['Max pixel value in x,y,z= ' num2str(max(newpoints(:,1))) ', ' num2str(max(newpoints(:,2))) ', ' num2str(max(newpoints(:,3)))]) %debug output
%disp(['Min pixel value in x,y,z= ' num2str(min(newpoints(:,1))) ', ' num2str(min(newpoints(:,2))) ', ' num2str(min(newpoints(:,3)))]) %debug output


filaments=cell(nfil,1);%a cell array in which each element contans an array with the coordinates of all points in a filament

table_out=[];


for i=1:nfil
    indices=point_IDs(i,point_IDs(i,:)~=-1)+1;%the indices of the points in this filament -- NOTE that the first point in the amira output has index 0, thus adding 1 here
    npoints=size(indices,2);%number of points in this filament
    filaments{i}=zeros(npoints,3);%initialising the array of coordinates for the points in this filament
    for j=1:npoints
       filaments{i}(j,:)=newpoints(indices(j),:);
    end
end

%disp(['There will be ' num2str(nfil) ' filaments'])
for i=1:nfil
    %i
    currfil=template_filament;
    currfil.radius=5;
    currfil.subunits_dphi=crop_dphi;
    currfil.subunits_dz=cropdist;
    currfil.backbone_interval=2;
    currfil.name=['mfilamentWithTorsion_' num2str(i)];
    currfil.file=[model_path 'mfilamentWithTorsion_' num2str(i) '.omd'];
    currfil.points=filaments{i};
    currfil.group_labels=ones(1,size(filaments{i},1));
    currfil.crop_points=[]; %if these are not empty, backboneUpdate will not fill them
    currfil.crop_angles=[]; %if these are not empty, backboneUpdate will not fill them
    currfil.backbone=[]; %if these are not empty, backboneUpdate will not fill them
    currfil.control_points=[]; %if these are not empty, backboneUpdate will not fill them 
    currfil.individual_labels=1:size(filaments{i},1);
    currfil.backboneUpdate(); %this is quivalent to the "smoothen backbone" step in the GUI
    currtable=currfil.grepTable();
    currtable(:,21)=i+offset_reg_nr;
    table_out=dtmerge({table_out,currtable},'linear_tags', 1);
    currfil.saveInCatalogue();
end

table_out(:,20)=tomo_nr;

