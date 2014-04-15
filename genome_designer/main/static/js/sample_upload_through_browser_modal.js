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
    this.renderSamplesAwaitingUpload();
    this.renderTemplateUploader();
    this.renderSamplesUploader();
  },

  renderSamplesAwaitingUpload: function() {
    $('#gd-samples-upload-through-browser-awaiting-upload').empty();

    var requestData = {
        projectUid: this.model.get('uid')
    };

    $.get('/_/samples/get_samples_awaiting_upload', requestData,
        function(response) {
          _.each(response.sampleFilenameList, function(sampleLabel) {
            $('#gd-samples-upload-through-browser-awaiting-upload').append(
                sampleLabel + ', ');
          });
        });
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
      autoUpload: true,
      formData: {projectUid: this.model.get('uid')}
    }).on('fileuploadadd', function (e, data) {
      this.clearSampleDataUploadEror();
      data.context = $('<div/>').appendTo('#files');
      $.each(data.files, function (index, file) {
        var node = $('<p/>')
            .append($('<span/>').text(file.name));
        if (!index) {
          node
              .append('<br>')
              .append(uploadButton.clone(true).data(data));
        }
        node.appendTo(data.context);
      });
    }).on('fileuploadprogressall', function (e, data) {
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
});
