# -*- coding: utf-8 -*-
############################################################
#
# Autogenerated by the KBase type compiler -
# any changes made here will be overwritten
#
############################################################

from __future__ import print_function
# the following is a hack to get the baseclient to import whether we're in a
# package or not. This makes pep8 unhappy hence the annotations.
try:
    # baseclient and this client are in a package
    from .baseclient import BaseClient as _BaseClient  # @UnusedImport
except:
    # no they aren't
    from baseclient import BaseClient as _BaseClient  # @Reimport
import time


class WsLargeDataIO(object):

    def __init__(
            self, url=None, timeout=30 * 60, user_id=None,
            password=None, token=None, ignore_authrc=False,
            trust_all_ssl_certificates=False,
            auth_svc='https://ci.kbase.us/services/auth/api/legacy/KBase/Sessions/Login',
            service_ver='release',
            async_job_check_time_ms=100, async_job_check_time_scale_percent=150, 
            async_job_check_max_time_ms=300000):
        if url is None:
            raise ValueError('A url is required')
        self._service_ver = service_ver
        self._client = _BaseClient(
            url, timeout=timeout, user_id=user_id, password=password,
            token=token, ignore_authrc=ignore_authrc,
            trust_all_ssl_certificates=trust_all_ssl_certificates,
            auth_svc=auth_svc,
            async_job_check_time_ms=async_job_check_time_ms,
            async_job_check_time_scale_percent=async_job_check_time_scale_percent,
            async_job_check_max_time_ms=async_job_check_max_time_ms)

    def _check_job(self, job_id):
        return self._client._check_job('WsLargeDataIO', job_id)

    def _save_objects_submit(self, params, context=None):
        return self._client._submit_job(
             'WsLargeDataIO.save_objects', [params],
             self._service_ver, context)

    def save_objects(self, params, context=None):
        """
        Save objects to the workspace. Saving over a deleted object undeletes
        it.
        :param params: instance of type "SaveObjectsParams" (Input parameters
           for the "save_objects" function. Required parameters: id - the
           numerical ID of the workspace. workspace - optional workspace name
           alternative to id. objects - the objects to save. The object
           provenance is automatically pulled from the SDK runner.) ->
           structure: parameter "id" of Long, parameter "workspace" of
           String, parameter "objects" of list of type "ObjectSaveData" (An
           object and associated data required for saving. Required
           parameters: type - the workspace type string for the object. Omit
           the version information to use the latest version. data_json_file
           - the path to a JSON file containing the object data. Optional
           parameters: One of an object name or id. If no name or id is
           provided the name will be set to 'auto' with the object id
           appended as a string, possibly with -\d+ appended if that object
           id already exists as a name. name - the name of the object. objid
           - the id of the object to save over. meta - arbitrary
           user-supplied metadata for the object, not to exceed 16kb; if the
           object type specifies automatic metadata extraction with the 'meta
           ws' annotation, and your metadata name conflicts, then your
           metadata will be silently overwritten. hidden - true if this
           object should not be listed when listing workspace objects.) ->
           structure: parameter "type" of String, parameter "data_json_file"
           of String, parameter "name" of String, parameter "objid" of Long,
           parameter "meta" of mapping from String to String, parameter
           "hidden" of type "boolean" (A boolean - 0 for false, 1 for true.
           @range (0, 1))
        :returns: instance of list of type "object_info" (Information about
           an object, including user provided metadata. objid - the numerical
           id of the object. name - the name of the object. type - the type
           of the object. save_date - the save date of the object. ver - the
           version of the object. saved_by - the user that saved or copied
           the object. wsid - the id of the workspace containing the object.
           workspace - the name of the workspace containing the object. chsum
           - the md5 checksum of the object. size - the size of the object in
           bytes. meta - arbitrary user-supplied metadata about the object.)
           -> tuple of size 11: parameter "objid" of Long, parameter "name"
           of String, parameter "type" of String, parameter "save_date" of
           String, parameter "version" of Long, parameter "saved_by" of
           String, parameter "wsid" of Long, parameter "workspace" of String,
           parameter "chsum" of String, parameter "size" of Long, parameter
           "meta" of mapping from String to String
        """
        job_id = self._save_objects_submit(params, context)
        async_job_check_time = self._client.async_job_check_time
        while True:
            time.sleep(async_job_check_time)
            async_job_check_time = (async_job_check_time *
                self._client.async_job_check_time_scale_percent / 100.0)
            if async_job_check_time > self._client.async_job_check_max_time:
                async_job_check_time = self._client.async_job_check_max_time
            job_state = self._check_job(job_id)
            if job_state['finished']:
                return job_state['result'][0]

    def _get_objects_submit(self, params, context=None):
        return self._client._submit_job(
             'WsLargeDataIO.get_objects', [params],
             self._service_ver, context)

    def get_objects(self, params, context=None):
        """
        Get objects from the workspace.
        :param params: instance of type "GetObjectsParams" (Input parameters
           for the "get_objects" function. Required parameters: object_refs -
           a list of object references in the form X/Y/Z, where X is the
           workspace name or id, Y is the object name or id, and Z is the
           (optional) object version. In general, always use ids rather than
           names if possible to avoid race conditions. Optional parameters:
           ignore_errors - ignore any errors that occur when fetching an
           object and instead insert a null into the returned list.) ->
           structure: parameter "objects" of list of type
           "ObjectSpecification" (An Object Specification (OS). Inherits from
           ObjectIdentity. Specifies which object, and which parts of that
           object, to retrieve from the Workspace Service. The fields wsid,
           workspace, objid, name, ver, and ref are identical to the
           ObjectIdentity fields. REFERENCE FOLLOWING: Reference following
           guarantees that a user that has access to an object can always see
           a) objects that are referenced inside the object and b) objects
           that are referenced in the object's provenance. This ensures that
           the user has visibility into the entire provenance of the object
           and the object's object dependencies (e.g. references). The user
           must have at least read access to the object specified in this SO,
           but need not have access to any further objects in the reference
           chain, and those objects may be deleted. Optional reference
           following fields: list<string> obj_ref_path - a path to the
           desired object from the object specified in this OS. In other
           words, the object specified in this OS is assumed to be accessible
           to the user, and the objects in the object path represent a chain
           of references to the desired object at the end of the object path.
           If the references are all valid, the desired object will be
           returned. OBJECT SUBSETS: When selecting a subset of an array in
           an object, the returned array is compressed to the size of the
           subset, but the ordering of the array is maintained. For example,
           if the array stored at the 'feature' key of a Genome object has
           4000 entries, and the object paths provided are: /feature/7
           /feature/3015 /feature/700 The returned feature array will be of
           length three and the entries will consist, in order, of the 7th,
           700th, and 3015th entries of the original array. Optional object
           subset fields: list<string> included - the portions of the object
           to include in the object subset. boolean strict_maps - if true,
           throw an exception if the subset specification traverses a
           non-existant map key (default false) boolean strict_arrays - if
           true, throw an exception if the subset specification exceeds the
           size of an array (default true)) -> structure: parameter
           "workspace" of String, parameter "wsid" of Long, parameter "name"
           of String, parameter "objid" of Long, parameter "ver" of Long,
           parameter "ref" of String, parameter "obj_ref_path" of list of
           String, parameter "included" of list of String, parameter
           "strict_maps" of type "boolean" (A boolean - 0 for false, 1 for
           true. @range (0, 1)), parameter "strict_arrays" of type "boolean"
           (A boolean - 0 for false, 1 for true. @range (0, 1)), parameter
           "ignore_errors" of type "boolean" (A boolean - 0 for false, 1 for
           true. @range (0, 1))
        :returns: instance of type "GetObjectsResults" (Results from the
           get_objects function. list<ObjectData> data - the returned
           objects.) -> structure: parameter "data" of list of type
           "ObjectData" (The data and supplemental info for an object.
           UnspecifiedObject data - the object's data or subset data.
           object_info info - information about the object.) -> structure:
           parameter "data_json_file" of String, parameter "info" of type
           "object_info" (Information about an object, including user
           provided metadata. objid - the numerical id of the object. name -
           the name of the object. type - the type of the object. save_date -
           the save date of the object. ver - the version of the object.
           saved_by - the user that saved or copied the object. wsid - the id
           of the workspace containing the object. workspace - the name of
           the workspace containing the object. chsum - the md5 checksum of
           the object. size - the size of the object in bytes. meta -
           arbitrary user-supplied metadata about the object.) -> tuple of
           size 11: parameter "objid" of Long, parameter "name" of String,
           parameter "type" of String, parameter "save_date" of String,
           parameter "version" of Long, parameter "saved_by" of String,
           parameter "wsid" of Long, parameter "workspace" of String,
           parameter "chsum" of String, parameter "size" of Long, parameter
           "meta" of mapping from String to String
        """
        job_id = self._get_objects_submit(params, context)
        async_job_check_time = self._client.async_job_check_time
        while True:
            time.sleep(async_job_check_time)
            async_job_check_time = (async_job_check_time *
                self._client.async_job_check_time_scale_percent / 100.0)
            if async_job_check_time > self._client.async_job_check_max_time:
                async_job_check_time = self._client.async_job_check_max_time
            job_state = self._check_job(job_id)
            if job_state['finished']:
                return job_state['result'][0]

    def status(self, context=None):
        job_id = self._client._submit_job('WsLargeDataIO.status', 
            [], self._service_ver, context)
        async_job_check_time = self._client.async_job_check_time
        while True:
            time.sleep(async_job_check_time)
            async_job_check_time = (async_job_check_time *
                self._client.async_job_check_time_scale_percent / 100.0)
            if async_job_check_time > self._client.async_job_check_max_time:
                async_job_check_time = self._client.async_job_check_max_time
            job_state = self._check_job(job_id)
            if job_state['finished']:
                return job_state['result'][0]