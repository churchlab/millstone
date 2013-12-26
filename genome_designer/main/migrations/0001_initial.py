# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'UserProfile'
        db.create_table(u'main_userprofile', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('user', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['auth.User'], unique=True)),
            ('uid', self.gf('django.db.models.fields.CharField')(default='3280f41b', unique=True, max_length=36)),
        ))
        db.send_create_signal(u'main', ['UserProfile'])

        # Adding model 'Dataset'
        db.create_table(u'main_dataset', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('uid', self.gf('django.db.models.fields.CharField')(default='4ccdf7e9', unique=True, max_length=36)),
            ('type', self.gf('django.db.models.fields.CharField')(max_length=40)),
            ('label', self.gf('django.db.models.fields.CharField')(max_length=256)),
            ('filesystem_location', self.gf('django.db.models.fields.CharField')(max_length=512, blank=True)),
            ('status', self.gf('django.db.models.fields.CharField')(default='READY', max_length=40)),
        ))
        db.send_create_signal(u'main', ['Dataset'])

        # Adding model 'Project'
        db.create_table(u'main_project', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('uid', self.gf('django.db.models.fields.CharField')(default='b7d66ef2', unique=True, max_length=36)),
            ('owner', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['main.UserProfile'])),
            ('title', self.gf('django.db.models.fields.CharField')(max_length=256)),
            ('s3_backed', self.gf('django.db.models.fields.BooleanField')(default=False)),
        ))
        db.send_create_signal(u'main', ['Project'])

        # Adding model 'ReferenceGenome'
        db.create_table(u'main_referencegenome', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('uid', self.gf('django.db.models.fields.CharField')(default='9f1473d5', unique=True, max_length=36)),
            ('project', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['main.Project'])),
            ('label', self.gf('django.db.models.fields.CharField')(max_length=256)),
            ('num_chromosomes', self.gf('django.db.models.fields.IntegerField')()),
            ('num_bases', self.gf('django.db.models.fields.BigIntegerField')()),
            ('variant_key_map', self.gf('jsonfield.fields.JSONField')()),
        ))
        db.send_create_signal(u'main', ['ReferenceGenome'])

        # Adding M2M table for field dataset_set on 'ReferenceGenome'
        m2m_table_name = db.shorten_name(u'main_referencegenome_dataset_set')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('referencegenome', models.ForeignKey(orm[u'main.referencegenome'], null=False)),
            ('dataset', models.ForeignKey(orm[u'main.dataset'], null=False))
        ))
        db.create_unique(m2m_table_name, ['referencegenome_id', 'dataset_id'])

        # Adding model 'ExperimentSample'
        db.create_table(u'main_experimentsample', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('uid', self.gf('django.db.models.fields.CharField')(default='3def9b65', unique=True, max_length=36)),
            ('project', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['main.Project'])),
            ('label', self.gf('django.db.models.fields.CharField')(max_length=256)),
            ('group', self.gf('django.db.models.fields.CharField')(max_length=256)),
            ('well', self.gf('django.db.models.fields.CharField')(max_length=256)),
            ('num_reads', self.gf('django.db.models.fields.BigIntegerField')(default=-1)),
        ))
        db.send_create_signal(u'main', ['ExperimentSample'])

        # Adding M2M table for field dataset_set on 'ExperimentSample'
        m2m_table_name = db.shorten_name(u'main_experimentsample_dataset_set')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('experimentsample', models.ForeignKey(orm[u'main.experimentsample'], null=False)),
            ('dataset', models.ForeignKey(orm[u'main.dataset'], null=False))
        ))
        db.create_unique(m2m_table_name, ['experimentsample_id', 'dataset_id'])

        # Adding model 'AlignmentGroup'
        db.create_table(u'main_alignmentgroup', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('uid', self.gf('django.db.models.fields.CharField')(default='8f5a90c1', unique=True, max_length=36)),
            ('label', self.gf('django.db.models.fields.CharField')(max_length=256, blank=True)),
            ('reference_genome', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['main.ReferenceGenome'])),
            ('aligner', self.gf('django.db.models.fields.CharField')(max_length=10)),
            ('start_time', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, blank=True)),
            ('end_time', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, blank=True)),
        ))
        db.send_create_signal(u'main', ['AlignmentGroup'])

        # Adding M2M table for field dataset_set on 'AlignmentGroup'
        m2m_table_name = db.shorten_name(u'main_alignmentgroup_dataset_set')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('alignmentgroup', models.ForeignKey(orm[u'main.alignmentgroup'], null=False)),
            ('dataset', models.ForeignKey(orm[u'main.dataset'], null=False))
        ))
        db.create_unique(m2m_table_name, ['alignmentgroup_id', 'dataset_id'])

        # Adding model 'ExperimentSampleToAlignment'
        db.create_table(u'main_experimentsampletoalignment', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('uid', self.gf('django.db.models.fields.CharField')(default='4d3c7226', unique=True, max_length=36)),
            ('alignment_group', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['main.AlignmentGroup'])),
            ('experiment_sample', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['main.ExperimentSample'])),
        ))
        db.send_create_signal(u'main', ['ExperimentSampleToAlignment'])

        # Adding M2M table for field dataset_set on 'ExperimentSampleToAlignment'
        m2m_table_name = db.shorten_name(u'main_experimentsampletoalignment_dataset_set')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('experimentsampletoalignment', models.ForeignKey(orm[u'main.experimentsampletoalignment'], null=False)),
            ('dataset', models.ForeignKey(orm[u'main.dataset'], null=False))
        ))
        db.create_unique(m2m_table_name, ['experimentsampletoalignment_id', 'dataset_id'])

        # Adding model 'Variant'
        db.create_table(u'main_variant', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('uid', self.gf('django.db.models.fields.CharField')(default='088ec42c', unique=True, max_length=36)),
            ('type', self.gf('django.db.models.fields.CharField')(max_length=40)),
            ('reference_genome', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['main.ReferenceGenome'])),
            ('chromosome', self.gf('django.db.models.fields.CharField')(max_length=256, blank=True)),
            ('position', self.gf('django.db.models.fields.BigIntegerField')()),
            ('ref_value', self.gf('django.db.models.fields.TextField')()),
        ))
        db.send_create_signal(u'main', ['Variant'])

        # Adding model 'VariantCallerCommonData'
        db.create_table(u'main_variantcallercommondata', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('variant', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['main.Variant'])),
            ('source_dataset', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['main.Dataset'])),
            ('data', self.gf('jsonfield.fields.JSONField')()),
        ))
        db.send_create_signal(u'main', ['VariantCallerCommonData'])

        # Adding model 'VariantAlternate'
        db.create_table(u'main_variantalternate', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('uid', self.gf('django.db.models.fields.CharField')(default='39c78a1f', max_length=36)),
            ('variant', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['main.Variant'], null=True)),
            ('alt_value', self.gf('django.db.models.fields.TextField')()),
            ('is_primary', self.gf('django.db.models.fields.BooleanField')(default=True)),
            ('data', self.gf('jsonfield.fields.JSONField')()),
        ))
        db.send_create_signal(u'main', ['VariantAlternate'])

        # Adding model 'VariantEvidence'
        db.create_table(u'main_variantevidence', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('uid', self.gf('django.db.models.fields.CharField')(default='d82922ce', max_length=36)),
            ('experiment_sample', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['main.ExperimentSample'])),
            ('variant_caller_common_data', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['main.VariantCallerCommonData'])),
            ('data', self.gf('jsonfield.fields.JSONField')()),
        ))
        db.send_create_signal(u'main', ['VariantEvidence'])

        # Adding M2M table for field variantalternate_set on 'VariantEvidence'
        m2m_table_name = db.shorten_name(u'main_variantevidence_variantalternate_set')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('variantevidence', models.ForeignKey(orm[u'main.variantevidence'], null=False)),
            ('variantalternate', models.ForeignKey(orm[u'main.variantalternate'], null=False))
        ))
        db.create_unique(m2m_table_name, ['variantevidence_id', 'variantalternate_id'])

        # Adding model 'VariantToVariantSet'
        db.create_table(u'main_varianttovariantset', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('variant', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['main.Variant'])),
            ('variant_set', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['main.VariantSet'])),
        ))
        db.send_create_signal(u'main', ['VariantToVariantSet'])

        # Adding M2M table for field sample_variant_set_association on 'VariantToVariantSet'
        m2m_table_name = db.shorten_name(u'main_varianttovariantset_sample_variant_set_association')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('varianttovariantset', models.ForeignKey(orm[u'main.varianttovariantset'], null=False)),
            ('experimentsample', models.ForeignKey(orm[u'main.experimentsample'], null=False))
        ))
        db.create_unique(m2m_table_name, ['varianttovariantset_id', 'experimentsample_id'])

        # Adding model 'VariantSet'
        db.create_table(u'main_variantset', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('uid', self.gf('django.db.models.fields.CharField')(default='da8ef027', unique=True, max_length=36)),
            ('label', self.gf('django.db.models.fields.CharField')(max_length=256)),
            ('reference_genome', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['main.ReferenceGenome'])),
        ))
        db.send_create_signal(u'main', ['VariantSet'])

        # Adding M2M table for field dataset_set on 'VariantSet'
        m2m_table_name = db.shorten_name(u'main_variantset_dataset_set')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('variantset', models.ForeignKey(orm[u'main.variantset'], null=False)),
            ('dataset', models.ForeignKey(orm[u'main.dataset'], null=False))
        ))
        db.create_unique(m2m_table_name, ['variantset_id', 'dataset_id'])

        # Adding model 'VariantFilter'
        db.create_table(u'main_variantfilter', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
        ))
        db.send_create_signal(u'main', ['VariantFilter'])

        # Adding model 'Region'
        db.create_table(u'main_region', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('uid', self.gf('django.db.models.fields.CharField')(default='f063e81e', unique=True, max_length=36)),
            ('reference_genome', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['main.ReferenceGenome'])),
            ('label', self.gf('django.db.models.fields.CharField')(max_length=256)),
            ('type', self.gf('django.db.models.fields.CharField')(max_length=40)),
        ))
        db.send_create_signal(u'main', ['Region'])

        # Adding model 'RegionInterval'
        db.create_table(u'main_regioninterval', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('region', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['main.Region'])),
            ('start', self.gf('django.db.models.fields.BigIntegerField')()),
            ('end', self.gf('django.db.models.fields.BigIntegerField')()),
        ))
        db.send_create_signal(u'main', ['RegionInterval'])

        # Adding model 'S3File'
        db.create_table(u'main_s3file', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('bucket', self.gf('django.db.models.fields.CharField')(max_length=200)),
            ('key', self.gf('django.db.models.fields.CharField')(max_length=200)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=200, null=True)),
            ('created_at', self.gf('django.db.models.fields.DateTimeField')(auto_now_add=True, blank=True)),
        ))
        db.send_create_signal(u'main', ['S3File'])


    def backwards(self, orm):
        # Deleting model 'UserProfile'
        db.delete_table(u'main_userprofile')

        # Deleting model 'Dataset'
        db.delete_table(u'main_dataset')

        # Deleting model 'Project'
        db.delete_table(u'main_project')

        # Deleting model 'ReferenceGenome'
        db.delete_table(u'main_referencegenome')

        # Removing M2M table for field dataset_set on 'ReferenceGenome'
        db.delete_table(db.shorten_name(u'main_referencegenome_dataset_set'))

        # Deleting model 'ExperimentSample'
        db.delete_table(u'main_experimentsample')

        # Removing M2M table for field dataset_set on 'ExperimentSample'
        db.delete_table(db.shorten_name(u'main_experimentsample_dataset_set'))

        # Deleting model 'AlignmentGroup'
        db.delete_table(u'main_alignmentgroup')

        # Removing M2M table for field dataset_set on 'AlignmentGroup'
        db.delete_table(db.shorten_name(u'main_alignmentgroup_dataset_set'))

        # Deleting model 'ExperimentSampleToAlignment'
        db.delete_table(u'main_experimentsampletoalignment')

        # Removing M2M table for field dataset_set on 'ExperimentSampleToAlignment'
        db.delete_table(db.shorten_name(u'main_experimentsampletoalignment_dataset_set'))

        # Deleting model 'Variant'
        db.delete_table(u'main_variant')

        # Deleting model 'VariantCallerCommonData'
        db.delete_table(u'main_variantcallercommondata')

        # Deleting model 'VariantAlternate'
        db.delete_table(u'main_variantalternate')

        # Deleting model 'VariantEvidence'
        db.delete_table(u'main_variantevidence')

        # Removing M2M table for field variantalternate_set on 'VariantEvidence'
        db.delete_table(db.shorten_name(u'main_variantevidence_variantalternate_set'))

        # Deleting model 'VariantToVariantSet'
        db.delete_table(u'main_varianttovariantset')

        # Removing M2M table for field sample_variant_set_association on 'VariantToVariantSet'
        db.delete_table(db.shorten_name(u'main_varianttovariantset_sample_variant_set_association'))

        # Deleting model 'VariantSet'
        db.delete_table(u'main_variantset')

        # Removing M2M table for field dataset_set on 'VariantSet'
        db.delete_table(db.shorten_name(u'main_variantset_dataset_set'))

        # Deleting model 'VariantFilter'
        db.delete_table(u'main_variantfilter')

        # Deleting model 'Region'
        db.delete_table(u'main_region')

        # Deleting model 'RegionInterval'
        db.delete_table(u'main_regioninterval')

        # Deleting model 'S3File'
        db.delete_table(u'main_s3file')


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
            'uid': ('django.db.models.fields.CharField', [], {'default': "'643365c5'", 'unique': 'True', 'max_length': '36'})
        },
        u'main.dataset': {
            'Meta': {'object_name': 'Dataset'},
            'filesystem_location': ('django.db.models.fields.CharField', [], {'max_length': '512', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'label': ('django.db.models.fields.CharField', [], {'max_length': '256'}),
            'status': ('django.db.models.fields.CharField', [], {'default': "'READY'", 'max_length': '40'}),
            'type': ('django.db.models.fields.CharField', [], {'max_length': '40'}),
            'uid': ('django.db.models.fields.CharField', [], {'default': "'a541ca0e'", 'unique': 'True', 'max_length': '36'})
        },
        u'main.experimentsample': {
            'Meta': {'object_name': 'ExperimentSample'},
            'dataset_set': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['main.Dataset']", 'null': 'True', 'blank': 'True'}),
            'group': ('django.db.models.fields.CharField', [], {'max_length': '256'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'label': ('django.db.models.fields.CharField', [], {'max_length': '256'}),
            'num_reads': ('django.db.models.fields.BigIntegerField', [], {'default': '-1'}),
            'project': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.Project']"}),
            'uid': ('django.db.models.fields.CharField', [], {'default': "'8a1075e0'", 'unique': 'True', 'max_length': '36'}),
            'well': ('django.db.models.fields.CharField', [], {'max_length': '256'})
        },
        u'main.experimentsampletoalignment': {
            'Meta': {'object_name': 'ExperimentSampleToAlignment'},
            'alignment_group': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.AlignmentGroup']"}),
            'dataset_set': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['main.Dataset']", 'null': 'True', 'blank': 'True'}),
            'experiment_sample': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.ExperimentSample']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'uid': ('django.db.models.fields.CharField', [], {'default': "'f6095264'", 'unique': 'True', 'max_length': '36'})
        },
        u'main.project': {
            'Meta': {'object_name': 'Project'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'owner': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.UserProfile']"}),
            's3_backed': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'title': ('django.db.models.fields.CharField', [], {'max_length': '256'}),
            'uid': ('django.db.models.fields.CharField', [], {'default': "'427ded09'", 'unique': 'True', 'max_length': '36'})
        },
        u'main.referencegenome': {
            'Meta': {'object_name': 'ReferenceGenome'},
            'dataset_set': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['main.Dataset']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'label': ('django.db.models.fields.CharField', [], {'max_length': '256'}),
            'num_bases': ('django.db.models.fields.BigIntegerField', [], {}),
            'num_chromosomes': ('django.db.models.fields.IntegerField', [], {}),
            'project': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.Project']"}),
            'uid': ('django.db.models.fields.CharField', [], {'default': "'c2aa1768'", 'unique': 'True', 'max_length': '36'}),
            'variant_key_map': ('jsonfield.fields.JSONField', [], {})
        },
        u'main.region': {
            'Meta': {'object_name': 'Region'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'label': ('django.db.models.fields.CharField', [], {'max_length': '256'}),
            'reference_genome': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.ReferenceGenome']"}),
            'type': ('django.db.models.fields.CharField', [], {'max_length': '40'}),
            'uid': ('django.db.models.fields.CharField', [], {'default': "'01598bab'", 'unique': 'True', 'max_length': '36'})
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
            'uid': ('django.db.models.fields.CharField', [], {'default': "'0122f58a'", 'unique': 'True', 'max_length': '36'}),
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
            'uid': ('django.db.models.fields.CharField', [], {'default': "'8a6a1cdc'", 'unique': 'True', 'max_length': '36'})
        },
        u'main.variantalternate': {
            'Meta': {'object_name': 'VariantAlternate'},
            'alt_value': ('django.db.models.fields.TextField', [], {}),
            'data': ('jsonfield.fields.JSONField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'is_primary': ('django.db.models.fields.BooleanField', [], {'default': 'True'}),
            'uid': ('django.db.models.fields.CharField', [], {'default': "'0ec6f115'", 'max_length': '36'}),
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
            'uid': ('django.db.models.fields.CharField', [], {'default': "'1afb584f'", 'max_length': '36'}),
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
            'uid': ('django.db.models.fields.CharField', [], {'default': "'4b2d9f60'", 'unique': 'True', 'max_length': '36'}),
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