# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):

        # Changing field 'Region.uid'
        db.alter_column(u'main_region', 'uid', self.gf('django.db.models.fields.CharField')(unique=True, max_length=8))

        # Changing field 'VariantAlternate.uid'
        db.alter_column(u'main_variantalternate', 'uid', self.gf('django.db.models.fields.CharField')(unique=True, max_length=8))
        # Adding unique constraint on 'VariantAlternate', fields ['uid']
        db.create_unique(u'main_variantalternate', ['uid'])


        # Changing field 'UserProfile.uid'
        db.alter_column(u'main_userprofile', 'uid', self.gf('django.db.models.fields.CharField')(unique=True, max_length=8))

        # Changing field 'VariantEvidence.uid'
        db.alter_column(u'main_variantevidence', 'uid', self.gf('django.db.models.fields.CharField')(unique=True, max_length=8))
        # Adding unique constraint on 'VariantEvidence', fields ['uid']
        db.create_unique(u'main_variantevidence', ['uid'])


        # Changing field 'Variant.uid'
        db.alter_column(u'main_variant', 'uid', self.gf('django.db.models.fields.CharField')(unique=True, max_length=8))

        # Changing field 'ExperimentSampleToAlignment.uid'
        db.alter_column(u'main_experimentsampletoalignment', 'uid', self.gf('django.db.models.fields.CharField')(unique=True, max_length=8))

        # Changing field 'VariantSet.uid'
        db.alter_column(u'main_variantset', 'uid', self.gf('django.db.models.fields.CharField')(unique=True, max_length=8))

        # Changing field 'AlignmentGroup.uid'
        db.alter_column(u'main_alignmentgroup', 'uid', self.gf('django.db.models.fields.CharField')(unique=True, max_length=8))

        # Changing field 'ReferenceGenome.uid'
        db.alter_column(u'main_referencegenome', 'uid', self.gf('django.db.models.fields.CharField')(unique=True, max_length=8))
        # Adding field 'Dataset.filesystem_idx_location'
        db.add_column(u'main_dataset', 'filesystem_idx_location',
                      self.gf('django.db.models.fields.CharField')(max_length=512, null=True, blank=True),
                      keep_default=False)


        # Changing field 'Dataset.uid'
        db.alter_column(u'main_dataset', 'uid', self.gf('django.db.models.fields.CharField')(unique=True, max_length=8))

        # Changing field 'ExperimentSample.uid'
        db.alter_column(u'main_experimentsample', 'uid', self.gf('django.db.models.fields.CharField')(unique=True, max_length=8))

        # Changing field 'Project.uid'
        db.alter_column(u'main_project', 'uid', self.gf('django.db.models.fields.CharField')(unique=True, max_length=8))

    def backwards(self, orm):
        # Removing unique constraint on 'VariantEvidence', fields ['uid']
        db.delete_unique(u'main_variantevidence', ['uid'])

        # Removing unique constraint on 'VariantAlternate', fields ['uid']
        db.delete_unique(u'main_variantalternate', ['uid'])


        # Changing field 'Region.uid'
        db.alter_column(u'main_region', 'uid', self.gf('django.db.models.fields.CharField')(max_length=36, unique=True))

        # Changing field 'VariantAlternate.uid'
        db.alter_column(u'main_variantalternate', 'uid', self.gf('django.db.models.fields.CharField')(max_length=36))

        # Changing field 'UserProfile.uid'
        db.alter_column(u'main_userprofile', 'uid', self.gf('django.db.models.fields.CharField')(max_length=36, unique=True))

        # Changing field 'VariantEvidence.uid'
        db.alter_column(u'main_variantevidence', 'uid', self.gf('django.db.models.fields.CharField')(max_length=36))

        # Changing field 'Variant.uid'
        db.alter_column(u'main_variant', 'uid', self.gf('django.db.models.fields.CharField')(max_length=36, unique=True))

        # Changing field 'ExperimentSampleToAlignment.uid'
        db.alter_column(u'main_experimentsampletoalignment', 'uid', self.gf('django.db.models.fields.CharField')(max_length=36, unique=True))

        # Changing field 'VariantSet.uid'
        db.alter_column(u'main_variantset', 'uid', self.gf('django.db.models.fields.CharField')(max_length=36, unique=True))

        # Changing field 'AlignmentGroup.uid'
        db.alter_column(u'main_alignmentgroup', 'uid', self.gf('django.db.models.fields.CharField')(max_length=36, unique=True))

        # Changing field 'ReferenceGenome.uid'
        db.alter_column(u'main_referencegenome', 'uid', self.gf('django.db.models.fields.CharField')(max_length=36, unique=True))
        # Deleting field 'Dataset.filesystem_idx_location'
        db.delete_column(u'main_dataset', 'filesystem_idx_location')


        # Changing field 'Dataset.uid'
        db.alter_column(u'main_dataset', 'uid', self.gf('django.db.models.fields.CharField')(max_length=36, unique=True))

        # Changing field 'ExperimentSample.uid'
        db.alter_column(u'main_experimentsample', 'uid', self.gf('django.db.models.fields.CharField')(max_length=36, unique=True))

        # Changing field 'Project.uid'
        db.alter_column(u'main_project', 'uid', self.gf('django.db.models.fields.CharField')(max_length=36, unique=True))

    models = {
        u'auth.group': {
            'Meta': {'object_name': 'Group'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '80'}),
            'permissions': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['auth.Permission']", 'symmetrical': 'False', 'blank': 'True'})
        },
        u'auth.permission': {
            'Meta': {'ordering': "(u'content_type__app_label', u'content_type__model', u'codename')", 'unique_together': "((u'content_type', u'codename'),)", 'object_name': 'Permission'},
            'codename': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'content_type': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['contenttypes.ContentType']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        u'auth.user': {
            'Meta': {'object_name': 'User'},
            'date_joined': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now'}),
            'email': ('django.db.models.fields.EmailField', [], {'max_length': '75', 'blank': 'True'}),
            'first_name': ('django.db.models.fields.CharField', [], {'max_length': '30', 'blank': 'True'}),
            'groups': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['auth.Group']", 'symmetrical': 'False', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'is_active': ('django.db.models.fields.BooleanField', [], {'default': 'True'}),
            'is_staff': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'is_superuser': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'last_login': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now'}),
            'last_name': ('django.db.models.fields.CharField', [], {'max_length': '30', 'blank': 'True'}),
            'password': ('django.db.models.fields.CharField', [], {'max_length': '128'}),
            'user_permissions': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['auth.Permission']", 'symmetrical': 'False', 'blank': 'True'}),
            'username': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '30'})
        },
        u'contenttypes.contenttype': {
            'Meta': {'ordering': "('name',)", 'unique_together': "(('app_label', 'model'),)", 'object_name': 'ContentType', 'db_table': "'django_content_type'"},
            'app_label': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'model': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '100'})
        },
        u'main.alignmentgroup': {
            'Meta': {'object_name': 'AlignmentGroup'},
            'aligner': ('django.db.models.fields.CharField', [], {'max_length': '10'}),
            'dataset_set': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['main.Dataset']", 'null': 'True', 'blank': 'True'}),
            'end_time': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'label': ('django.db.models.fields.CharField', [], {'max_length': '256', 'blank': 'True'}),
            'reference_genome': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.ReferenceGenome']"}),
            'start_time': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'}),
            'uid': ('django.db.models.fields.CharField', [], {'default': "'0288b3e5'", 'unique': 'True', 'max_length': '8'})
        },
        u'main.dataset': {
            'Meta': {'object_name': 'Dataset'},
            'filesystem_idx_location': ('django.db.models.fields.CharField', [], {'max_length': '512', 'null': 'True', 'blank': 'True'}),
            'filesystem_location': ('django.db.models.fields.CharField', [], {'max_length': '512', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'label': ('django.db.models.fields.CharField', [], {'max_length': '256'}),
            'status': ('django.db.models.fields.CharField', [], {'default': "'READY'", 'max_length': '40'}),
            'type': ('django.db.models.fields.CharField', [], {'max_length': '40'}),
            'uid': ('django.db.models.fields.CharField', [], {'default': "'6676ee28'", 'unique': 'True', 'max_length': '8'})
        },
        u'main.experimentsample': {
            'Meta': {'object_name': 'ExperimentSample'},
            'dataset_set': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['main.Dataset']", 'null': 'True', 'blank': 'True'}),
            'group': ('django.db.models.fields.CharField', [], {'max_length': '256'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'label': ('django.db.models.fields.CharField', [], {'max_length': '256'}),
            'num_reads': ('django.db.models.fields.BigIntegerField', [], {'default': '-1'}),
            'project': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.Project']"}),
            'uid': ('django.db.models.fields.CharField', [], {'default': "'275fb523'", 'unique': 'True', 'max_length': '8'}),
            'well': ('django.db.models.fields.CharField', [], {'max_length': '256'})
        },
        u'main.experimentsampletoalignment': {
            'Meta': {'object_name': 'ExperimentSampleToAlignment'},
            'alignment_group': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.AlignmentGroup']"}),
            'dataset_set': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['main.Dataset']", 'null': 'True', 'blank': 'True'}),
            'experiment_sample': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.ExperimentSample']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'uid': ('django.db.models.fields.CharField', [], {'default': "'575f061e'", 'unique': 'True', 'max_length': '8'})
        },
        u'main.project': {
            'Meta': {'object_name': 'Project'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'owner': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.UserProfile']"}),
            's3_backed': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'title': ('django.db.models.fields.CharField', [], {'max_length': '256'}),
            'uid': ('django.db.models.fields.CharField', [], {'default': "'82bc43fe'", 'unique': 'True', 'max_length': '8'})
        },
        u'main.referencegenome': {
            'Meta': {'object_name': 'ReferenceGenome'},
            'dataset_set': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['main.Dataset']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'label': ('django.db.models.fields.CharField', [], {'max_length': '256'}),
            'num_bases': ('django.db.models.fields.BigIntegerField', [], {}),
            'num_chromosomes': ('django.db.models.fields.IntegerField', [], {}),
            'project': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.Project']"}),
            'uid': ('django.db.models.fields.CharField', [], {'default': "'7c17c23b'", 'unique': 'True', 'max_length': '8'}),
            'variant_key_map': ('jsonfield.fields.JSONField', [], {})
        },
        u'main.region': {
            'Meta': {'object_name': 'Region'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'label': ('django.db.models.fields.CharField', [], {'max_length': '256'}),
            'reference_genome': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.ReferenceGenome']"}),
            'type': ('django.db.models.fields.CharField', [], {'max_length': '40'}),
            'uid': ('django.db.models.fields.CharField', [], {'default': "'f43dbed9'", 'unique': 'True', 'max_length': '8'})
        },
        u'main.regioninterval': {
            'Meta': {'object_name': 'RegionInterval'},
            'end': ('django.db.models.fields.BigIntegerField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'region': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.Region']"}),
            'start': ('django.db.models.fields.BigIntegerField', [], {})
        },
        u'main.s3file': {
            'Meta': {'object_name': 'S3File'},
            'bucket': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'created_at': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'key': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200', 'null': 'True'})
        },
        u'main.userprofile': {
            'Meta': {'object_name': 'UserProfile'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'uid': ('django.db.models.fields.CharField', [], {'default': "'e0aae88d'", 'unique': 'True', 'max_length': '8'}),
            'user': ('django.db.models.fields.related.OneToOneField', [], {'to': u"orm['auth.User']", 'unique': 'True'})
        },
        u'main.variant': {
            'Meta': {'object_name': 'Variant'},
            'chromosome': ('django.db.models.fields.CharField', [], {'max_length': '256', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'position': ('django.db.models.fields.BigIntegerField', [], {}),
            'ref_value': ('django.db.models.fields.TextField', [], {}),
            'reference_genome': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.ReferenceGenome']"}),
            'type': ('django.db.models.fields.CharField', [], {'max_length': '40'}),
            'uid': ('django.db.models.fields.CharField', [], {'default': "'fc29b81c'", 'unique': 'True', 'max_length': '8'})
        },
        u'main.variantalternate': {
            'Meta': {'object_name': 'VariantAlternate'},
            'alt_value': ('django.db.models.fields.TextField', [], {}),
            'data': ('jsonfield.fields.JSONField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'is_primary': ('django.db.models.fields.BooleanField', [], {'default': 'True'}),
            'uid': ('django.db.models.fields.CharField', [], {'default': "'6e421b28'", 'unique': 'True', 'max_length': '8'}),
            'variant': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.Variant']", 'null': 'True'})
        },
        u'main.variantcallercommondata': {
            'Meta': {'object_name': 'VariantCallerCommonData'},
            'data': ('jsonfield.fields.JSONField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'source_dataset': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.Dataset']"}),
            'variant': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.Variant']"})
        },
        u'main.variantevidence': {
            'Meta': {'object_name': 'VariantEvidence'},
            'data': ('jsonfield.fields.JSONField', [], {}),
            'experiment_sample': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.ExperimentSample']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'uid': ('django.db.models.fields.CharField', [], {'default': "'35a4d172'", 'unique': 'True', 'max_length': '8'}),
            'variant_caller_common_data': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.VariantCallerCommonData']"}),
            'variantalternate_set': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['main.VariantAlternate']", 'symmetrical': 'False'})
        },
        u'main.variantfilter': {
            'Meta': {'object_name': 'VariantFilter'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'})
        },
        u'main.variantset': {
            'Meta': {'object_name': 'VariantSet'},
            'dataset_set': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['main.Dataset']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'label': ('django.db.models.fields.CharField', [], {'max_length': '256'}),
            'reference_genome': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.ReferenceGenome']"}),
            'uid': ('django.db.models.fields.CharField', [], {'default': "'3eb9d8ad'", 'unique': 'True', 'max_length': '8'}),
            'variants': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['main.Variant']", 'null': 'True', 'through': u"orm['main.VariantToVariantSet']", 'blank': 'True'})
        },
        u'main.varianttovariantset': {
            'Meta': {'object_name': 'VariantToVariantSet'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'sample_variant_set_association': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['main.ExperimentSample']", 'null': 'True', 'blank': 'True'}),
            'variant': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.Variant']"}),
            'variant_set': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.VariantSet']"})
        }
    }

    complete_apps = ['main']