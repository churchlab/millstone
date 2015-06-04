/**
 * Component for modal that facilitates uploading samples through the browser.
 *
 * Usage of jQuery File Upload based on:
 * https://github.com/sigurdga/django-jquery-file-upload
 */


gd.SampleUploadThroughBrowserModal = Backbone.View.extend({

  initialize: function() {
    if (!this.model) {
      throw "gd.SampleUploadThroughBrowserModal requires model.";
    }

    this.render();
  },

  render: function() {
    // Register model change listener
    this.model.on('change:samplesAwaitingUpload', _.bind(function() {
      $('#gd-samples-upload-through-browser-awaiting-upload').empty();
      _.each(this.model.get('samplesAwaitingUpload'),
          _.bind(function(sampleLabel) {
            $('#gd-samples-upload-through-browser-awaiting-upload').append(
                sampleLabel + ', ');
          }, this));
      }, this));

    this.renderSamplesAwaitingUpload();
    this.renderTemplateUploader();
    this.renderSamplesUploader();
  },

  renderSamplesAwaitingUpload: function() {
    var requestData = {
        projectUid: this.model.get('uid')
    };

    $.get('/_/samples/get_samples_awaiting_upload', requestData,
        _.bind(function(response) {
          this.model.set({
              'samplesAwaitingUpload': response.sampleFilenameList
          });
        }, this));
  },

  renderTemplateUploader: function() {
    // Create the uploader component.
    this.templateUploader = $('#gd-samples-template-fileupload').fileupload({
      url: '/_/samples/samples_upload_through_browser_template',
      autoUpload: true,
      formData: {projectUid: this.model.get('uid')}
    }).on('fileuploadadd', _.bind(function (e, data) {
      this.clearTemplateUploadError();
    }, this)).on('fileuploaddone', _.bind(function(e, data) {
      var result = data.response().result;
      if ('error' in result) {
        this.showTemplateUploadError('ERROR: ' + result.error);
      } else {
        this.renderSamplesAwaitingUpload();
      }
    }, this));

  },

  renderSamplesUploader: function() {
    // Upload button added when a file is added for upload.
    var uploadButton = $('<button/>')
        .addClass('btn btn-primary')
        .prop('disabled', true)
        .text('Processing...')
        .on('click', function () {
          var $this = $(this),
          data = $this.data();
          $this
              .off('click')
              .text('Abort')
              .on('click', function () {
                $this.remove();
                data.abort();
              });
          data.submit().always(function () {
            $this.remove();
          });
        });

    // Create the uploader component.
    this.samplesUploader = $('#gd-samples-fastq-fileupload').fileupload({
      url: '/_/samples/samples_upload_through_browser_sample_data',
      autoUpload: false,
      formData: {projectUid: this.model.get('uid')}
    }).on('fileuploadadd', _.bind(function (e, data) {
      this.clearSampleDataUploadEror();
      data.context = $('<div/>').appendTo('#files');
      for (var index = 0; index < data.files.length; index++) {
        var file = data.files[index];

        // Check that this file is expected.
        if (!_.contains(this.model.get('samplesAwaitingUpload'), file.name)) {
          this.showSampleDataUploadError('Unexpected file ' + file.name +
              '. Please make sure you registered it correctly using the ' +
              'template above.');
          continue;
        }

        // Append to ui.
        var node = $('<p/>')
            .append($('<span/>').text(file.name));
        if (!index) {
          node
              .append('<br>')
              .append(uploadButton.clone(true).data(data));
        }
        node.appendTo(data.context);

        // Kick off upload.
        data.process().done(function () {
          data.submit();
        });
      }
    }, this)).on('fileuploadprogressall', function (e, data) {
      var progress = parseInt(data.loaded / data.total * 100, 10);
      $('#progress .progress-bar').css(
        'width',
        progress + '%'
        );
    }).on('fileuploaddone', _.bind(function (e, data) {
      // Maybe show an error.
      var result = data.response().result;
      if ('error' in result) {
        this.showSampleDataUploadError(result.error);
      } else {
        this.clearSampleDataUploadEror();
      }

      // Update the child node buttons.
      $.each(data.files, function (index, file) {
        var statusEl = $(data.context.children()[index]);
        var loadingButton = statusEl.children('button');
        loadingButton.text('Done');
      });

      // Update the files awaiting upload.
      this.renderSamplesAwaitingUpload();
    }, this)).on('fileuploadfail', function (e, data) {
      $.each(data.files, function (index, file) {
        var error = $('<span/>').text(file.error);
        $(data.context.children()[index])
            .append('<br>')
            .append(error);
      });
    });
  },

  /**
   * Shows an error when there is a problem with uploading a sample template.
   */
  showTemplateUploadError: function(errorMsg) {
    this.clearTemplateUploadError();
    var errorEl = $('#gd-samples-upload-through-browser-modal-template-error');
    errorEl.append('<p>' + errorMsg + '</p');
  },

  /** Clears the error message. */
  clearTemplateUploadError: function() {
    var errorEl = $('#gd-samples-upload-through-browser-modal-template-error');
    errorEl.empty();
  },

  /**
   * Shows an error when there is a problem with uploading a sample.
   *
   * TODO: Make this show up next to the correct file status div.
   */
  showSampleDataUploadError: function(errorMsg) {
    this.clearSampleDataUploadEror();
    var errorEl = $('#gd-samples-upload-through-browser-modal-error');
    errorEl.append('<p>' + errorMsg + '</p');
  },

  /** Clears the error message. */
  clearSampleDataUploadEror: function() {
    var errorEl = $('#gd-samples-upload-through-browser-modal-error');
    errorEl.empty();
  },

  destroy: function() {
    this.remove();
    this.unbind();
  }
});
