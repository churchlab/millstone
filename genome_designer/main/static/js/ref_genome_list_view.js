/**
 * @fileoverview Reference Genome List view.
 */


gd.RefGenomeListView = Backbone.View.extend({
  el: '#gd-page-container',

  initialize: function() {
    this.render();
    this.decorate_new_button();
  },

  render: function() {

    $('#gd-sidenav-link-refgenomes').addClass('active');

    this.datatable = new gd.DataTableComponent({
        el: $('#gd-ref-genome-list-view-datatable-hook'),
        objList: REF_GENOME_LIST_DATA['obj_list'],
        fieldConfig: REF_GENOME_LIST_DATA['field_config']
    });

    this.uploader = this.$("#uploadDiv").fineUploaderS3({
      debug: true,
      request: {
        endpoint: this.$("#uploadDiv").data("endpoint"),
        accessKey: this.$("#uploadDiv").data("accesskey")
      },
      signature: {
        endpoint: this.$("#uploadDiv").data("signature")
      },
      uploadSuccess: {
        endpoint: this.$("#uploadDiv").data("success")
      },
      objectProperties: {
        key: $.proxy(function(fileId) {
          var filename = this.uploader.fineUploader("getName", fileId);
          return "uploads/" + qq.getUniqueId() + "-" + filename;
        }, this)
      },
      retry: {
        enableAuto: true
      },
      chunking: {
        enabled: true
      },
      deleteFile: {
        endpoint: this.$("#uploadDiv").data("delete"),
        enabled: true,
        forceConfirm: true
      },
      callbacks: {
        onError: function(id, name, reason) {
          alert(reason);
        }
      }
    }).on('complete', $.proxy(function(id, name, response, xhr) {
      var sid = xhr.s3file_id;
      $.post(this.$("#uploadDiv").data("import"), 
            {
              's3file_id': sid,
              'refGenomeLabel': this.$("#uploadDiv").find("#refGenomeLabel").val(),
              'importFileFormat': $("#uploadDiv").find("input[name=importFileFormat]:checked").val(),
            },
            function(data) {
            }
      );}, this));
  },

  decorate_new_button: function() {

    html_string = (
      '<div class="btn-group">' +
      '  <a class="btn dropdown-toggle btn-primary" data-toggle="dropdown" href="#">' +
      '    New' +
      '    <span class="caret"></span>' +
      '  </a>' +
      '  <ul class="dropdown-menu">' +
      '    <li role="presentation"><a role="menuitem" tabindex="-1" href="#modalFromFile" data-toggle="modal">From Server Location...</a></li>'
    );

    if (IS_S3_BACKEND) {
      html_string = html_string + 
        '    <li role="presentation"><a role="menuitem" tabindex="-1" href="#modalUpload" data-toggle="modal">From S3 Bucket...</a></li>';
    };
    html_string = html_string + 
      '    <li role="presentation"><a role="menuitem" tabindex="-1" href="#modalFromNCBI" data-toggle="modal">From NCBI Accession...</a></li>' +
      '  </ul>' +
      '</div>'

    $("div.gd-new-button").html(html_string);
  }

});
