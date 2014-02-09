@CeleryManager = @CeleryManager || {}
@CeleryManager.config = @CeleryManager.config || {}
@CeleryManager.classes = @CeleryManager.classes || {}

class CeleryManager.HomeView extends Backbone.View
  initialize: =>
    @template = _.template $("#home-view-template").text()
    @render()

  success: (header, message) =>
    @$("#message").text message
    @$("#message-header").text header
    @$("#message-box").removeClass("alert-danger hidden").addClass("alert-success")
    @$("#save-button").attr "disabled", null

  error: (header, message) =>
    @$("#message").text message
    @$("#message-header").text header
    @$("#message-box").removeClass("alert-success hidden").addClass("alert-danger").show()
    @$("#save-button").attr "disabled", disabled

  render: =>
    @$el.append @template()

    jQuery.validator.addMethod("domainOrIpv4", (value, element, param) ->
        ipv4regex = /^(25[0-5]|2[0-4]\d|[01]?\d\d?)\.(25[0-5]|2[0-4]\d|[01]?\d\d?)\.(25[0-5]|2[0-4]\d|[01]?\d\d?)\.(25[0-5]|2[0-4]\d|[01]?\d\d?)$/i
        ec2regex = /^ec2-(25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)-(25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)-(25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)-(25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)\.compute-1\.amazonaws\.com$/i
        return this.optional(element) || ec2regex.test(value) || ipv4regex.test(value)
      , "Please enter a valid IPv4 or domain name of EC2 instance.")
    @$("form").validate rules:
      host: 
        required: true,
        domainOrIpv4: true
      user: 
        required: true,
      password: 
        required: false

    @$("#test-button").click =>
      $.ajax
        url: CeleryManager.config.urlRoot + 'test'
        type: 'post'
        data: @$("form").serialize()
        success: (data, textStatus, jqXHR) =>
          if data.status == "error"
            @error "Failed the test!", data.message
          else
            @success "Passed the test!", "Click Save to continue"
      return false

    @$("#save-button").click =>
      $.ajax
        url: CeleryManager.config.urlRoot + 'save'
        type: 'post'
        data: @$("form").serialize()
        success: (data, textStatus, jqXHR) =>
          if data.status == "success"
            @success "Saved!", ""
      return false

class CeleryManager.Router extends Backbone.Router
  routes:
    ""  :   "home"

  initialize: =>
    jQuery.validator.setDefaults
      errorPlacement: (error, element) ->
        if(element.parent().hasClass('input-prepend') || element.parent().hasClass('input-append'))
          error.insertAfter(element.parent())
        else
          error.insertAfter(element)
      errorElement: "small"
      wrapper: "div"
      highlight: (element) ->
        $(element).closest('.form-group').addClass('has-error')
      success: (element) ->
        $(element).closest('.form-group').removeClass('has-error')

    CeleryManager.config.urlRoot = window.location.href.split("?")[0]

  home: =>
    @currentView = new CeleryManager.HomeView el: $("#startup-view")

$ ->
  CeleryManager.app = new CeleryManager.Router()
  Backbone.history.start()
