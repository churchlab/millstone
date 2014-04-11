/**
 * Component for modal that facilitates uploading samples through the browser.
 */


gd.SampleUploadThroughBrowserModal = Backbone.View.extend({

  initialize: function() {
    if (!this.model) {
      throw "gd.SampleUploadThroughBrowserModal requires model.";
    }

    this.render();
  },

  render: function() {
    this.renderTemplateUploader();
    this.renderSamplesUploader();
  },

  renderTemplateUploader: function() {
    if (this.templateUploader) {
      return;
    }

    // Create the uploader component.
    this.templateUploader = $('#gd-samples-template-fileupload').fileupload({
      url: '/_/samples/samples_upload_through_browser_template',
      autoUpload: true,
      formData: {projectUid: this.model.get('uid')}
    });

    // Handle the result.
    this.templateUploader.on('fileuploaddone', function (e, data) {
      var result = data.response().result;
      if ('error' in result) {
        console.log(result.error);
      } else {
        window.location.reload();
      }
    });

  },

  renderSamplesUploader: function() {
    if (this.samplesUploader) {
      return;
    }

    // Create the uploader component.
    this.samplesUploader = $('#gd-samples-fastq-fileupload').fileupload({
      url: '/_/samples/samples_upload_through_browser_sample_data',
      autoUpload: true,
      formData: {projectUid: this.model.get('uid')}
    });

    // Upload button.
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

    // Provide ui affordance for file added to upload queue.
    this.samplesUploader.on('fileuploadadd', function (e, data) {
      data.context = $('<div/>').appendTo('#files');
      $.each(data.files, function (index, file) {
          var node = $('<p/>').append($('<span/>').text(file.name));
          if (!index) {
              node.append('<br>')
                  .append(uploadButton.clone(true).data(data));
          }
          node.appendTo(data.context);
      });
    })

    // Update progress bar.
    this.samplesUploader.on('fileuploadprogressall', function (e, data) {
        var progress = parseInt(data.loaded / data.total * 100, 10);
        $('#progress .progress-bar').css(
            'width',
            progress + '%'
        );
    });

    // Handle the response.
    this.samplesUploader.on('fileuploaddone', _.bind(function (e, data) {
      var result = data.response().result;
      if ('error' in result) {
        this.showSampleDataUploadError(result.error);
      } else {
        this.clearSampleDataUploadEror();
      }
    }, this));
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
