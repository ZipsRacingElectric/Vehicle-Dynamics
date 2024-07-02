CM4SL and TM4SL models can be updated using the matlab scripts ModelUpdate_CM12.m
and ModelUpdate_TM12.m in your Simulink folder (src_cm4sl/src_tm4sl. The scripts
update all Simulink models in the folder and save them as<YourModel>_updated.mdl.
A report with all changes is saved as <YourModel>_changes.html. If an update
fails, the matlab error is shown in the matlab command window and and the script
proceeds with the next model.

Updating may fail if you inserted models at unpredicted places or changed the
model structure. In that case you can try to comment out the line
update_<component_xyz> which is shown in the error mesasage and redo the update
for the specfic model. To restrict the update to a single model select the
matching model(s) for the variable FileList. The skipped part has then to be
updated manually. All files used for the update are stored in the folder
ModelUpdate so you can also try to adjust them to your models if errros occur.

As an alternative you can update your Simulink model or parts of it manually.
Copy the new version of the appropriate Simulink model from the Simulink folder
(src_cm4sl/src_tm4sl) in the project directory and copy your changes into this
new model. If you need support for upgrading your CM4SL/TM4SL models, do not
hesitate to contact the CarMaker Support Team (support@ipg-automotive.com).