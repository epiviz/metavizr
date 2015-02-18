/**
 * Created by: Florin Chelaru
 * Date: 10/2/13
 * Time: 12:46 PM
 */

goog.provide('epiviz.EpiViz');

goog.require('epiviz.Config');
goog.require('epiviz.data.DataProviderFactory');
goog.require('epiviz.data.DataManager');
goog.require('epiviz.data.DataProvider');


/**
 * @param {epiviz.Config} config
 * @param {epiviz.ui.LocationManager} locationManager
 * @param {epiviz.measurements.MeasurementsManager} measurementsManager
 * @param {epiviz.ui.ControlManager} controlManager
 * @param {epiviz.data.DataManager} dataManager
 * @param {epiviz.ui.charts.ChartFactory} chartFactory
 * @param {epiviz.ui.charts.ChartManager} chartManager
 * @param {epiviz.workspaces.WorkspaceManager} workspaceManager
 * @param {epiviz.workspaces.UserManager} userManager
 * @param {epiviz.ui.WebArgsManager} webArgsManager
 * @constructor
 */
epiviz.EpiViz = function(config, locationManager, measurementsManager, controlManager, dataManager, chartFactory, chartManager, workspaceManager, userManager, webArgsManager) {
  /**
   * @type {epiviz.Config}
   * @private
   */
  this._config = config;

  /**
   * @type {epiviz.ui.LocationManager}
   * @private
   */
  this._locationManager = locationManager;

  /**
   * @type {epiviz.measurements.MeasurementsManager}
   * @private
   */
  this._measurementsManager = measurementsManager;

  /**
   * @type {epiviz.ui.ControlManager}
   * @private
   */
  this._controlManager = controlManager;

  /**
   * @type {epiviz.data.DataManager}
   * @private
   */
  this._dataManager = dataManager;

  /**
   * @type {epiviz.ui.charts.ChartFactory}
   * @private
   */
  this._chartFactory = chartFactory;

  /**
   * @type {epiviz.ui.charts.ChartManager}
   * @private
   */
  this._chartManager = chartManager;

  /**
   * @type {epiviz.workspaces.WorkspaceManager}
   * @private
   */
  this._workspaceManager = workspaceManager;

  /**
   * @type {epiviz.workspaces.UserManager}
   * @private
   */
  this._userManager = userManager;

  /**
   * @type {epiviz.ui.WebArgsManager}
   * @private
   */
  this._webArgsManager = webArgsManager;

  // Register for UI events

  this._registerRequestSeqInfos();
  this._registerRequestMeasurements();

  this._registerUiAddChart();
  this._registerUiSaveWorkspace();
  this._registerUiDeleteActiveWorkspace();
  this._registerUiLoginLinkClicked();
  this._registerUiSearchWorkspaces();
  this._registerUiActiveWorkspaceChanged();
  this._registerUiSearch();

  // Register for Data events

  this._registerDataAddMeasurements();
  this._registerDataRemoveMeasurements();
  this._registerDataAddChart();
  this._registerDataRemoveChart();
  this._registerDataAddSeqInfos();
  this._registerDataRemoveSeqNames();
  this._registerDataNavigate();
  this._registerDataRedraw();
  this._registerDataGetCurrentLocation();

  // Register for Workspace events

  this._registerRequestWorkspaces();
  this._registerWorkspacesLoaded();
  this._registerActiveWorkspaceChanged();
  this._registerLocationChanged();

  /*
   * Prevent closing if workspace has changed
   */
  var self = this;
  window.onbeforeunload = function() {
    if (epiviz.workspaces.UserManager.USER_STATUS.loggedIn && self._workspaceManager.activeWorkspace().changed()) {
      return 'There are unsaved changes in the current workspace. Do you wish to discard them?';
    }
    return undefined;
   };
};

/**
 * A map of settings that are used as input for the EpiViz configuration
 * @type {*}
 */
epiviz.EpiViz.SETTINGS = {};

/**
 * @type {string}
 * @const
 */
epiviz.EpiViz.VERSION = '2';

epiviz.EpiViz.prototype.start = function() {
  this._controlManager.initialize();

  this._workspaceManager.initialize();

  this._measurementsManager.initialize();

  this._locationManager.initialize();
};

/**
 * @returns {epiviz.Config}
 */
epiviz.EpiViz.prototype.config = function() {
  return this._config;
};

/**
 * @param {epiviz.ui.charts.ChartType} type
 * @param {epiviz.measurements.MeasurementSet} measurements
 * @param {string} [chartId] If specified, then this will be
 *   the id of the newly created chart. Otherwise, a new one
 *   will be generated.
 * @param {epiviz.ui.charts.ChartProperties} [chartProperties]
 * @returns {string} the id of the chart just created
 * @private
 */
epiviz.EpiViz.prototype._addChart = function(type, measurements, chartId, chartProperties) {
  chartId = this._chartManager.addChart(type, measurements, chartId, chartProperties);
  var range = this._workspaceManager.activeWorkspace().range();
  this._chartManager.addChartsLoaderAnimation(chartId);
  var chartMeasurementsMap = {};
  chartMeasurementsMap[chartId] = measurements;
  var self = this;
  this._dataManager.getData(range, chartMeasurementsMap,
    function(chartId, data) {
      self._chartManager.updateCharts(range, data, [chartId]);
    });

  return chartId;
};

/*****************************************************************************
 * UI Events                                                                 *
 *****************************************************************************/

/**
 * @private
 */
epiviz.EpiViz.prototype._registerRequestSeqInfos = function() {
  var self = this;
  this._locationManager.onRequestSeqInfos().addListener(new epiviz.events.EventListener(
    function() {
      self._dataManager.getSeqInfos(function(seqInfos) {
        self._locationManager.updateSeqInfos(seqInfos);
      });
    }));
};

/**
 * @private
 */
epiviz.EpiViz.prototype._registerRequestMeasurements = function() {
  var self = this;
  this._measurementsManager.onRequestMeasurements().addListener(new epiviz.events.EventListener(
    function() {
      self._dataManager.getMeasurements(function(measurements) {
        self._measurementsManager.addMeasurements(measurements);

        //self._workspaceManager.initialize();
      });
    }));
};

/**
 * @private
 */
epiviz.EpiViz.prototype._registerRequestWorkspaces = function() {
  var self = this;
  this._workspaceManager.onRequestWorkspaces().addListener(new epiviz.events.EventListener(
    /**
     * @param {{activeWorkspaceId: string}} e
     */
    function(e) {
      self._dataManager.getWorkspaces(function(rawWorkspaces) {
        var ws = [];
        var activeWorkspace = null;
        for (var i = 0; i < rawWorkspaces.length; ++i) {
           var w = epiviz.workspaces.Workspace.fromRawObject(rawWorkspaces[i], self._chartFactory);

          if (w.id() === null) {
            // This is a workspace retrieved using e.activeWorkspaceId
            // and belonging to another user
            activeWorkspace = w;
            continue;
          }

          if (w.id() == e.activeWorkspaceId) {
            activeWorkspace = w;
          }

          ws.push(w);
        }

        self._workspaceManager.updateWorkspaces(ws, activeWorkspace, e.activeWorkspaceId);
        self._workspaceManager.activeWorkspace().resetChanged();
       }, '', e.activeWorkspaceId);
    }));
};

/**
 * @private
 */
epiviz.EpiViz.prototype._registerUiAddChart = function() {
  var self = this;
  this._controlManager.onAddChart().addListener(new epiviz.events.EventListener(
    /** @param {{type: epiviz.ui.charts.ChartType, measurements: epiviz.measurements.MeasurementSet}} e */
    function(e) {
      self._addChart(e.type, e.measurements);
    }));
};

/**
 * @private
 */
epiviz.EpiViz.prototype._registerUiSaveWorkspace = function() {
  var self = this;
  this._controlManager.onSaveWorkspace().addListener(new epiviz.events.EventListener(
    /**
     * @param {{name: string, id: ?string}} e
     */
    function(e) {
      // If there is a workspace with this name in the set of workspaces, overwrite it
      // Otherwise, create a copy of the active workspace with the new name, and save that.
      var workspace = self._workspaceManager.getByName(e.name);
      if (workspace) {
        workspace = self._workspaceManager.activeWorkspace().copy(e.name, workspace.id());
      } else {
        workspace = self._workspaceManager.activeWorkspace().copy(e.name);
      }

      self._dataManager.saveWorkspace(workspace, function(id) {
        workspace = workspace.copy(workspace.name(), id);
        workspace.resetChanged();

        self._workspaceManager.updateWorkspace(workspace);
        self._workspaceManager.changeActiveWorkspace(id);
      });
    }
  ));
};

/**
 * @private
 */
epiviz.EpiViz.prototype._registerUiDeleteActiveWorkspace = function() {
  var self = this;
  this._controlManager.onDeleteActiveWorkspace().addListener(new epiviz.events.EventListener(
    function() {
      self._dataManager.deleteWorkspace(self._workspaceManager.activeWorkspace());
      self._workspaceManager.deleteActiveWorkspace();
    }
  ));
};

/**
 * @private
 */
epiviz.EpiViz.prototype._registerUiLoginLinkClicked = function() {
  var self = this;
  this._controlManager.onLoginLinkClicked().addListener(new epiviz.events.EventListener(function() {
    self._userManager.toggleLogin();
  }));
};

/**
 * @private
 */
epiviz.EpiViz.prototype._registerUiSearchWorkspaces = function() {
  var self = this;
  this._controlManager.onSearchWorkspaces().addListener(new epiviz.events.EventListener(
    /**
     * @param {{searchTerm: string, callback: function(Array.<epiviz.workspaces.Workspace>)}} e
     */
    function(e) {
      self._dataManager.getWorkspaces(function(workspaces) {
        e.callback(workspaces);
      }, e.searchTerm, e.searchTerm);
    }
  ))
};

/**
 * @private
 */
epiviz.EpiViz.prototype._registerUiActiveWorkspaceChanged = function() {
  var self = this;
  this._controlManager.onActiveWorkspaceChanged().addListener(new epiviz.events.EventListener(
    /**
     * @param {{oldValue: {id: string, name: string}, newValue: {id: string, name: string}, cancel: function}} e
     */
    function(e) {

      var doChangeActiveWorkspace = function() {
        if (!self._workspaceManager.get(e.newValue.id)) {
          // The requested workspace id belongs to another user, so it has to be retrieved
          self._dataManager.getWorkspaces(function(rawWorkspaces) {
            var result = null;
            for (var i = 0; i < rawWorkspaces.length; ++i) {
              var w = epiviz.workspaces.Workspace.fromRawObject(rawWorkspaces[i], self._chartFactory);

              if (w.id() === null) {
                // This is a workspace retrieved using e.activeWorkspaceId
                // and belonging to another user
                result = w;
                break;
              }
            }

            if (result) {
              self._workspaceManager.changeActiveWorkspace(e.newValue.id, result);
            }
          }, e.newValue.name, e.newValue.id);
        } else { self._workspaceManager.changeActiveWorkspace(e.newValue.id); }
      };

      if (epiviz.workspaces.UserManager.USER_STATUS.loggedIn && !self._workspaceManager.activeWorkspaceChanging() && self._workspaceManager.activeWorkspace().changed()) {
        var dialog = new epiviz.ui.controls.MessageDialog(
          'Discard workspace changes',
          {
            Yes: function() {
              doChangeActiveWorkspace();
            },
            No: function() {
              e.cancel();
            }
          },
          'There are unsaved changes in the current workspace. Do you wish to discard them?',
          epiviz.ui.controls.MessageDialog.Icon.QUESTION);
        dialog.show();
      } else {
        doChangeActiveWorkspace();
      }
    }
  ));
};

/**
 * @private
 */
epiviz.EpiViz.prototype._registerUiSearch = function() {
  var self = this;
  this._controlManager.onSearch().addListener(new epiviz.events.EventListener(
    /**
     * @param {{searchTerm: string, callback: (function(Array.<{probe: string, gene: string, seqName: string, start: number, end: number}>))}} e
     */
    function(e) {
      self._dataManager.search(function(results) {
        e.callback(results);
      }, e.searchTerm);
    }));
};

/*****************************************************************************
 * Data                                                                      *
 *****************************************************************************/

/**
 * @private
 */
epiviz.EpiViz.prototype._registerDataAddMeasurements = function() {
  var self = this;
  this._dataManager.onRequestAddMeasurements().addListener(new epiviz.events.EventListener(
    /** @param {{measurements: epiviz.measurements.MeasurementSet, result: epiviz.events.EventResult}} e */
    function(e) {
      try {
        self._measurementsManager.addMeasurements(e.measurements);
        e.result.success = true;
      } catch (error) {
        e.result.success = false;
        e.result.errorMessage = error.toString();
      }
    }));
};

/**
 * @private
 */
epiviz.EpiViz.prototype._registerDataRemoveMeasurements = function() {
  var self = this;
  this._dataManager.onRequestRemoveMeasurements().addListener(new epiviz.events.EventListener(
    /** @param {{measurements: epiviz.measurements.MeasurementSet, result: epiviz.events.EventResult}} e */
      function(e) {
      try {
        self._measurementsManager.removeMeasurements(e.measurements);
        e.result.success = true;
      } catch (error) {
        e.result.success = false;
        e.result.errorMessage = error.toString();
      }
    }));
};

/**
 * @private
 */
epiviz.EpiViz.prototype._registerDataAddChart = function() {
  var self = this;
  this._dataManager.onRequestAddChart().addListener(new epiviz.events.EventListener(
    /** @param {{type: string, measurements: epiviz.measurements.MeasurementSet, result: epiviz.events.EventResult}} e */
    function(e) {
      try {
        var chartType = self._chartFactory.get(e.type);
        var chartId = self._addChart(chartType, e.measurements);
        e.result.success = true;
        e.result.value = { id: chartId };
      } catch (error) {
        e.result.success = false;
        e.result.errorMessage = error.toString();
      }
    }));
};

/**
 * @private
 */
epiviz.EpiViz.prototype._registerDataRemoveChart = function() {
  var self = this;
  this._dataManager.onRequestRemoveChart().addListener(new epiviz.events.EventListener(
    /**
     * @param {{id: string, result: epiviz.events.EventResult}} e
     */
    function(e) {
      try {
        self._chartManager.removeChart(e.id);
        e.result.success = true;
      } catch (error) {
        e.result.success = false;
        e.errorMessage = error.toString();
      }
    }));
};

/**
 * @private
 */
epiviz.EpiViz.prototype._registerDataAddSeqInfos = function() {
  var self = this;
  this._dataManager.onRequestAddSeqInfos().addListener(new epiviz.events.EventListener(
    /**
     * @param {{seqInfos: Array.<epiviz.datatypes.SeqInfo>, result: epiviz.events.EventResult}} e
     */
    function(e) {
      try {
        self._locationManager.addSeqInfos(e.seqInfos);
        e.result.success = true;
      } catch (error) {
        e.result.success = false;
        e.errorMessage = error.toString();
      }
    }));
};

/**
 * @private
 */
epiviz.EpiViz.prototype._registerDataRemoveSeqNames = function() {
  var self = this;
  this._dataManager.onRequestRemoveSeqNames().addListener(new epiviz.events.EventListener(
    /**
     * @param {{seqNames: Array.<string>, result: epiviz.events.EventResult}} e
     */
    function(e) {
      try {
        self._locationManager.removeSeqNames(e.seqNames);
        e.result.success = true;
      } catch (error) {
        e.result.success = false;
        e.errorMessage = error.toString();
      }
    }));
};

/**
 * @private
 */
epiviz.EpiViz.prototype._registerDataNavigate = function() {
  var self = this;
  this._dataManager.onRequestNavigate().addListener(new epiviz.events.EventListener(
    /**
     * @param {{range: epiviz.datatypes.GenomicRange, result: epiviz.events.EventResult}} e
     */
    function(e) {
      try {
        self._locationManager.changeCurrentLocation(e.range);
        e.result.success = true;
      } catch (error) {
        e.result.success = false;
        e.errorMessage = error.toString();
      }
    }));
};

/**
 * @private
 */
epiviz.EpiViz.prototype._registerDataRedraw = function() {
  var self = this;
  this._dataManager.onRequestRedraw().addListener(new epiviz.events.EventListener(
    /**
     * @param {{result: epiviz.events.EventResult}} e
     */
    function(e) {
      try {
        var currentLocation = self._locationManager.currentLocation();
        self._locationManager.changeCurrentLocation(currentLocation);
        e.result.success = true;
      } catch (error) {
        e.result.success = false;
        e.errorMessage = error.toString();
      }
    }));
};

/**
 * @private
 */
epiviz.EpiViz.prototype._registerDataGetCurrentLocation = function() {
  var self = this;
  this._dataManager.onRequestCurrentLocation().addListener(new epiviz.events.EventListener(
    /**
     * @param {{result: epiviz.events.EventResult}} e
     */
      function(e) {
      try {
        var currentLocation = self._locationManager.currentLocation();
        e.result.value = { seqName: currentLocation.seqName(), start: currentLocation.start(), end: currentLocation.end() };
        e.result.success = true;
      } catch (error) {
        e.result.success = false;
        e.errorMessage = error.toString();
      }
    }));
};

/*****************************************************************************
 * Workspaces                                                                *
 *****************************************************************************/

/**
 * @private
 */
epiviz.EpiViz.prototype._registerWorkspacesLoaded = function() {
  var self = this;
  this._workspaceManager.onWorkspacesLoaded().addListener(new epiviz.events.EventListener(
    /**
     * @param {{
     *   activeWorkspace: epiviz.workspaces.Workspace,
     *   workspaces: Array.<epiviz.workspaces.Workspace>
     * }} e
     */
    function(e) {
    }));
};

/**
 * @private
 */
epiviz.EpiViz.prototype._registerActiveWorkspaceChanged = function() {
  var self = this;
  this._workspaceManager.onActiveWorkspaceChanged().addListener(new epiviz.events.EventListener(
    /**
     * @param {{oldValue: epiviz.workspaces.Workspace, newValue: epiviz.workspaces.Workspace, workspaceId: string}} e
     */
    function(e) {
      self._workspaceManager.startChangingActiveWorkspace();

      self._controlManager.updateSelectedWorkspace({id: e.newValue.id(), name: e.newValue.name()});
      self._locationManager.changeCurrentLocation(e.newValue.range());

      self._measurementsManager.removeMeasurements(self._measurementsManager.computedMeasurements());
      self._measurementsManager.addMeasurements(e.newValue.computedMeasurements());

      self._chartManager.clear();

      /**
       * @type {Object.<epiviz.ui.charts.ChartType.DisplayType, Array.<{id: string, type: epiviz.ui.charts.ChartType, properties: epiviz.ui.charts.ChartProperties}>>}
       */
      var charts = e.newValue.charts();

      for (var displayType in charts) {
        if (!charts.hasOwnProperty(displayType)) { continue; }

        for (var i = 0; i < charts[displayType].length; ++i) {
          self._addChart(charts[displayType][i].type, charts[displayType][i].properties.measurements, charts[displayType][i].id, charts[displayType][i].properties.copy());
        }
      }

      self._workspaceManager.endChangingActiveWorkspace();
    }
  ));
};

/**
 * @private
 */
epiviz.EpiViz.prototype._registerLocationChanged = function() {
  var self = this;
  this._locationManager.onCurrentLocationChanged().addListener(new epiviz.events.EventListener(
    /**
     * @param {{oldValue: epiviz.datatypes.GenomicRange, newValue: epiviz.datatypes.GenomicRange}} e
     */
    function(e) {
      self._chartManager.addChartsLoaderAnimation();

      /** @type {Object.<string, epiviz.measurements.MeasurementSet>} */
      var chartMeasurementsMap = self._chartManager.chartsMeasurements();

      self._dataManager.getData(e.newValue, chartMeasurementsMap,
        function(chartId, data) {
          self._chartManager.updateCharts(e.newValue, data, [chartId]);
        });
    }));
};
