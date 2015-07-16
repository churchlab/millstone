***************
Troubleshooting
***************

-  I can log in via SSH but the web interface doesn't load!

    You've probably forgotten to allow access to your instance through
    web interfaces. This can be fixed by adding the following
    connections to your security group: \* All ICMP \* All TCP \* All
    UDP You can do this by going to the Network & Security -> Security
    Groups section of the EC2 dashboard and editing the security group
    that you created in your instance. If you've forgotten this can be
    found in the main instance dash on the far right under security
    groups. Click on that and you should be able to edit inbound rules
    by right clicking on the Group ID

-  I've managed to load the webpage but get a 502 bad gateway error!

    Millstone is probably loading up, try again in a few minutes.

-  Registration is closed.

    Only one user is allowed to register (as soon as the server boots
    up), and afterwards registration is closed.

-  Millstone just sits there after importing a template file.

    This could be any number of things. If your template file is
    formatted correctly, it could be a completely out of space error, so
    check that you've got room on your drive containing Millstone. File
    formatting is often the biggest problem in this stage, so be careful
    that you've escaped spaces in file names.

-  I want to make sure everything's going right, where can I find the
   logs?

    The logs are by default at /var/log/supervisor